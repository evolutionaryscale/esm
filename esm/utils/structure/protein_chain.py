from __future__ import annotations

import io
from dataclasses import asdict, dataclass, replace
from functools import cached_property
from pathlib import Path
from typing import Sequence, TypeVar, Union

import biotite.structure as bs
import brotli
import msgpack
import msgpack_numpy
import numpy as np
import torch
from Bio.Data import PDBData
from biotite.application.dssp import DsspApp
from biotite.database import rcsb
from biotite.structure.io.npz import NpzFile
from biotite.structure.io.pdb import PDBFile
from cloudpathlib import CloudPath
from scipy.spatial.distance import pdist, squareform
from torch import Tensor

from esm.utils import residue_constants as RC
from esm.utils.constants import esm3 as C
from esm.utils.misc import slice_python_object_as_numpy
from esm.utils.structure.affine3d import Affine3D
from esm.utils.structure.aligner import Aligner
from esm.utils.structure.metrics import compute_lddt_ca
from esm.utils.structure.normalize_coordinates import (
    apply_frame_to_coords,
    get_protein_normalization_frame,
    normalize_coordinates,
)

msgpack_numpy.patch()

CHAIN_ID_CONST = "A"


ArrayOrTensor = TypeVar("ArrayOrTensor", np.ndarray, Tensor)
PathLike = Union[str, Path, CloudPath]
PathOrBuffer = Union[PathLike, io.StringIO]


def index_by_atom_name(
    atom37: ArrayOrTensor, atom_names: str | list[str], dim: int = -2
) -> ArrayOrTensor:
    squeeze = False
    if isinstance(atom_names, str):
        atom_names = [atom_names]
        squeeze = True
    indices = [RC.atom_order[atom_name] for atom_name in atom_names]
    dim = dim % atom37.ndim
    index = tuple(slice(None) if dim != i else indices for i in range(atom37.ndim))
    result = atom37[index]  # type: ignore
    if squeeze:
        result = result.squeeze(dim)
    return result


def infer_CB(C, N, Ca, L: float = 1.522, A: float = 1.927, D: float = -2.143):
    """
    Inspired by a util in trDesign:
    https://github.com/gjoni/trDesign/blob/f2d5930b472e77bfacc2f437b3966e7a708a8d37/02-GD/utils.py#L92

    input:  3 coords (a,b,c), (L)ength, (A)ngle, and (D)ihedral
    output: 4th coord
    """
    norm = lambda x: x / np.sqrt(np.square(x).sum(-1, keepdims=True) + 1e-8)
    with np.errstate(invalid="ignore"):  # inf - inf = nan is ok here
        vec_bc = N - Ca
        vec_ba = N - C
    bc = norm(vec_bc)
    n = norm(np.cross(vec_ba, bc))
    m = [bc, np.cross(n, bc), n]
    d = [L * np.cos(A), L * np.sin(A) * np.cos(D), -L * np.sin(A) * np.sin(D)]
    return Ca + sum([m * d for m, d in zip(m, d)])


class AtomIndexer:
    def __init__(self, structure: ProteinChain, property: str, dim: int):
        self.structure = structure
        self.property = property
        self.dim = dim

    def __getitem__(self, atom_names: str | list[str]) -> np.ndarray:
        return index_by_atom_name(
            getattr(self.structure, self.property), atom_names, self.dim
        )


@dataclass
class ProteinChain:
    """Dataclass with atom37 representation of a single protein chain."""

    id: str
    sequence: str
    chain_id: str  # author chain id
    entity_id: int | None
    residue_index: np.ndarray
    insertion_code: np.ndarray
    atom37_positions: np.ndarray
    atom37_mask: np.ndarray
    confidence: np.ndarray

    def __post_init__(self):
        self.atom37_mask = self.atom37_mask.astype(bool)
        assert self.atom37_positions.shape[0] == len(self.sequence), (
            self.atom37_positions.shape,
            len(self.sequence),
        )
        assert self.atom37_mask.shape[0] == len(self.sequence), (
            self.atom37_mask.shape,
            len(self.sequence),
        )
        assert self.residue_index.shape[0] == len(self.sequence), (
            self.residue_index.shape,
            len(self.sequence),
        )
        assert self.insertion_code.shape[0] == len(self.sequence), (
            self.insertion_code.shape,
            len(self.sequence),
        )
        assert self.confidence.shape[0] == len(self.sequence), (
            self.confidence.shape,
            len(self.sequence),
        )

    @cached_property
    def atoms(self) -> AtomIndexer:
        return AtomIndexer(self, property="atom37_positions", dim=-2)

    @cached_property
    def atom_mask(self) -> AtomIndexer:
        return AtomIndexer(self, property="atom37_mask", dim=-1)

    @cached_property
    def atom_array(self) -> bs.AtomArray:
        atoms = []
        for res_name, res_idx, ins_code, positions, mask, conf in zip(
            self.sequence,
            self.residue_index,
            self.insertion_code,
            self.atom37_positions,
            self.atom37_mask.astype(bool),
            self.confidence,
        ):
            for i, pos in zip(np.where(mask)[0], positions[mask]):
                atom = bs.Atom(
                    coord=pos,
                    chain_id="A" if self.chain_id is None else self.chain_id,
                    res_id=res_idx,
                    ins_code=ins_code,
                    res_name=RC.restype_1to3.get(res_name, "UNK"),
                    hetero=False,
                    atom_name=RC.atom_types[i],
                    element=RC.atom_types[i][0],
                    b_factor=conf,
                )
                atoms.append(atom)
        return bs.array(atoms)

    @cached_property
    def residue_index_no_insertions(self) -> np.ndarray:
        return self.residue_index + np.cumsum(self.insertion_code != "")

    @cached_property
    def atom_array_no_insertions(self) -> bs.AtomArray:
        atoms = []
        for res_idx, (res_name, positions, mask, conf) in enumerate(
            zip(
                self.sequence,
                self.atom37_positions,
                self.atom37_mask.astype(bool),
                self.confidence,
            )
        ):
            for i, pos in zip(np.where(mask)[0], positions[mask]):
                atom = bs.Atom(
                    coord=pos,
                    # hard coded to as we currently only support single chain structures
                    chain_id=CHAIN_ID_CONST,
                    res_id=res_idx + 1,
                    res_name=RC.restype_1to3.get(res_name, "UNK"),
                    hetero=False,
                    atom_name=RC.atom_types[i],
                    element=RC.atom_types[i][0],
                    b_factor=conf,
                )
                atoms.append(atom)
        return bs.array(atoms)

    def __getitem__(self, idx: int | list[int] | slice | np.ndarray):
        if isinstance(idx, int):
            idx = [idx]

        sequence = slice_python_object_as_numpy(self.sequence, idx)
        return replace(
            self,
            sequence=sequence,
            residue_index=self.residue_index[..., idx],
            insertion_code=self.insertion_code[..., idx],
            atom37_positions=self.atom37_positions[..., idx, :, :],
            atom37_mask=self.atom37_mask[..., idx, :],
            confidence=self.confidence[..., idx],
        )

    def __len__(self):
        return len(self.sequence)

    def cbeta_contacts(self, distance_threshold: float = 8.0) -> np.ndarray:
        distance = self.pdist_CB
        contacts = (distance < distance_threshold).astype(np.int64)
        contacts[np.isnan(distance)] = -1
        np.fill_diagonal(contacts, -1)
        return contacts

    def to_npz(self, path: PathOrBuffer):
        f = NpzFile()
        f.set_structure(self.atom_array)
        f.write(path)

    def to_npz_string(self):
        f = NpzFile()
        f.set_structure(self.atom_array)
        buf = io.BytesIO()
        f.write(buf)
        return buf.getvalue()

    def to_structure_encoder_inputs(
        self,
        should_normalize_coordinates: bool = True,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        coords = torch.tensor(self.atom37_positions, dtype=torch.float32)
        plddt = torch.tensor(self.confidence, dtype=torch.float32)
        residue_index = torch.tensor(self.residue_index, dtype=torch.long)

        if should_normalize_coordinates:
            coords = normalize_coordinates(coords)
        return coords.unsqueeze(0), plddt.unsqueeze(0), residue_index.unsqueeze(0)

    def to_pdb(self, path: PathOrBuffer, include_insertions: bool = True):
        """Dssp works better w/o insertions."""
        f = PDBFile()
        if not include_insertions:
            f.set_structure(self.atom_array_no_insertions)
        else:
            f.set_structure(self.atom_array)
        f.write(path)

    def to_pdb_string(self, include_insertions: bool = True) -> str:
        buf = io.StringIO()
        self.to_pdb(buf, include_insertions=include_insertions)
        buf.seek(0)
        return buf.read()

    def state_dict(self, backbone_only=False):
        """This state dict is optimized for storage, so it turns things to fp16 whenever
        possible. Note that we also only support int32 residue indices, I'm hoping we don't
        need more than 2**32 residues..."""
        dct = {k: v for k, v in asdict(self).items()}
        for k, v in dct.items():
            if isinstance(v, np.ndarray):
                match v.dtype:
                    case np.int64:
                        dct[k] = v.astype(np.int32)
                    case np.float64 | np.float32:
                        dct[k] = v.astype(np.float16)
                    case _:
                        pass
        if backbone_only:
            dct["atom37_mask"][:, 3:] = False
        dct["atom37_positions"] = dct["atom37_positions"][dct["atom37_mask"]]
        return dct

    def to_blob(self, backbone_only=False) -> bytes:
        return brotli.compress(msgpack.dumps(self.state_dict(backbone_only)))

    @classmethod
    def from_state_dict(cls, dct):
        atom37 = np.full((*dct["atom37_mask"].shape, 3), np.nan)
        atom37[dct["atom37_mask"]] = dct["atom37_positions"]
        dct["atom37_positions"] = atom37
        dct = {
            k: (v.astype(np.float32) if k in ["atom37_positions", "confidence"] else v)
            for k, v in dct.items()
        }
        return cls(**dct)

    @classmethod
    def from_blob(cls, input: Path | str | io.BytesIO | bytes):
        """NOTE: blob + sparse coding + brotli + fp16 reduces memory
        of chains from 52G/1M chains to 20G/1M chains, I think this is a good first
        shot at compressing and dumping chains to disk. I'm sure there's better ways."""
        match input:
            case Path() | str():
                bytes = Path(input).read_bytes()
            case io.BytesIO():
                bytes = input.getvalue()
            case _:
                bytes = input
        return cls.from_state_dict(msgpack.loads(brotli.decompress(bytes)))

    def dssp(self):
        dssp = DsspApp.annotate_sse(self.atom_array_no_insertions)
        full_dssp = np.full(len(self.sequence), "X", dtype="<U1")
        full_dssp[self.atom37_mask.any(-1)] = dssp
        return full_dssp

    def sasa(self):
        arr = self.atom_array_no_insertions
        sasa_per_atom = bs.sasa(arr)  # type: ignore
        # Sum per-atom SASA into residue "bins", with np.bincount.
        assert arr.res_id is not None
        assert np.array_equal(
            np.sort(np.unique(arr.res_id)), np.arange(1, arr.res_id.max() + 1)
        ), "SASA calculation expected contiguous res_ids in range(1, len(chain)+1)"
        # NOTE: arr.res_id is 1-indexed, but np.bincount returns a sum for bin 0, so we strip.
        sasa_per_residue = np.bincount(arr.res_id, weights=sasa_per_atom)[1:]
        assert len(sasa_per_residue) == len(self)
        return sasa_per_residue

    def align(
        self,
        target: ProteinChain,
        mobile_inds: list[int] | np.ndarray | None = None,
        target_inds: list[int] | np.ndarray | None = None,
        only_use_backbone: bool = False,
    ):
        """
        Aligns the current protein to the provided target.

        Args:
            target (ProteinChain): The target protein to align to.
            mobile_inds (list[int], np.ndarray, optional): The indices of the mobile atoms to align. These are NOT residue indices
            target_inds (list[int], np.ndarray, optional): The indices of the target atoms to align. These are NOT residue indices
            only_use_backbone (bool, optional): If True, only align the backbone atoms.
        """
        aligner = Aligner(
            self if mobile_inds is None else self[mobile_inds],
            target if target_inds is None else target[target_inds],
            only_use_backbone,
        )

        return aligner.apply(self)

    def rmsd(
        self,
        target: ProteinChain,
        also_check_reflection: bool = False,
        mobile_inds: list[int] | np.ndarray | None = None,
        target_inds: list[int] | np.ndarray | None = None,
        only_compute_backbone_rmsd: bool = False,
    ):
        """
        Compute the RMSD between this protein chain and another.

        Args:
            target (ProteinChain): The target (other) protein chain to compare to.
            also_check_reflection (bool, optional): If True, also check if the reflection of the mobile atoms has a lower RMSD.
            mobile_inds (list[int], optional): The indices of the mobile atoms to align. These are NOT residue indices
            target_inds (list[int], optional): The indices of the target atoms to align. These are NOT residue indices
            only_compute_backbone_rmsd (bool, optional): If True, only compute the RMSD of the backbone atoms.
        """
        if isinstance(target, bs.AtomArray):
            raise ValueError(
                "Support for bs.AtomArray removed, use "
                "ProteinChain.from_atomarry for ProteinChain."
            )
        aligner = Aligner(
            self if mobile_inds is None else self[mobile_inds],
            target if target_inds is None else target[target_inds],
            only_compute_backbone_rmsd,
        )
        avg_rmsd = aligner.rmsd

        if not also_check_reflection:
            return avg_rmsd

        aligner = Aligner(
            self if mobile_inds is None else self[mobile_inds],
            target if target_inds is None else target[target_inds],
            only_compute_backbone_rmsd,
            use_reflection=True,
        )
        avg_rmsd_neg = aligner.rmsd

        return min(avg_rmsd, avg_rmsd_neg)

    def lddt_ca(
        self,
        native: ProteinChain,
        mobile_inds: list[int] | np.ndarray | None = None,
        target_inds: list[int] | np.ndarray | None = None,
        **kwargs,
    ) -> float | np.ndarray:
        """Compute the LDDT between this protein chain and another.
        NOTE: LDDT IS NOT SYMMETRIC. The call should always be prediction.lddt_ca(native).

        Arguments:
            native (ProteinChain): The ground truth protein chain
            mobile_inds (list[int], np.ndarray, optional): The indices of the mobile atoms to align. These are NOT residue indices
            target_inds (list[int], np.ndarray, optional): The indices of the target atoms to align. These are NOT residue indices

        Returns:
            float | np.ndarray: The LDDT score between the two protein chains, either
                a single float or per-residue LDDT scores if `per_residue` is True.
        """
        lddt = compute_lddt_ca(
            torch.tensor(self.atom37_positions[mobile_inds]).unsqueeze(0),
            torch.tensor(native.atom37_positions[target_inds]).unsqueeze(0),
            torch.tensor(native.atom37_mask[mobile_inds]).unsqueeze(0),
            **kwargs,
        )
        return float(lddt) if lddt.numel() == 1 else lddt.numpy().flatten()

    @classmethod
    def from_atom37(
        cls,
        atom37_positions: np.ndarray | torch.Tensor,
        *,
        id: str | None = None,
        sequence: str | None = None,
        chain_id: str | None = None,
        entity_id: int | None = None,
        residue_index: np.ndarray | torch.Tensor | None = None,
        insertion_code: np.ndarray | None = None,
        confidence: np.ndarray | torch.Tensor | None = None,
    ):
        if isinstance(atom37_positions, torch.Tensor):
            atom37_positions = atom37_positions.cpu().numpy()
            if atom37_positions.ndim == 4:
                if atom37_positions.shape[0] != 1:
                    raise ValueError(
                        f"Cannot handle batched inputs, atom37_positions has shape {atom37_positions.shape}"
                    )
                atom37_positions = atom37_positions[0]

        assert isinstance(atom37_positions, np.ndarray)
        seqlen = atom37_positions.shape[0]

        atom_mask = np.isfinite(atom37_positions).all(-1)

        if id is None:
            id = ""

        if sequence is None:
            sequence = "A" * seqlen

        if chain_id is None:
            chain_id = "A"

        if residue_index is None:
            residue_index = np.arange(1, seqlen + 1)
        elif isinstance(residue_index, torch.Tensor):
            residue_index = residue_index.cpu().numpy()
            assert isinstance(residue_index, np.ndarray)
            if residue_index.ndim == 2:
                if residue_index.shape[0] != 1:
                    raise ValueError(
                        f"Cannot handle batched inputs, residue_index has shape {residue_index.shape}"
                    )
                residue_index = residue_index[0]
        assert isinstance(residue_index, np.ndarray)

        if insertion_code is None:
            insertion_code = np.array(["" for _ in range(seqlen)])

        if confidence is None:
            confidence = np.ones(seqlen, dtype=np.float32)
        elif isinstance(confidence, torch.Tensor):
            confidence = confidence.cpu().numpy()
            assert isinstance(confidence, np.ndarray)
            if confidence.ndim == 2:
                if confidence.shape[0] != 1:
                    raise ValueError(
                        f"Cannot handle batched inputs, confidence has shape {confidence.shape}"
                    )
                confidence = confidence[0]
        assert isinstance(confidence, np.ndarray)

        return cls(
            id=id,
            sequence=sequence,
            chain_id=chain_id,
            entity_id=entity_id,
            atom37_positions=atom37_positions,
            atom37_mask=atom_mask,
            residue_index=residue_index,
            insertion_code=insertion_code,
            confidence=confidence,
        )

    @classmethod
    def from_backbone_atom_coordinates(
        cls,
        backbone_atom_coordinates: np.ndarray | torch.Tensor,
        **kwargs,
    ):
        """Create a ProteinChain from a set of backbone atom coordinates.

        This function simply expands the seqlen x 3 x 3 array of backbone atom
        coordinates to a seqlen x 37 x 3 array of all atom coordinates, with the padded
        positions set to infinity. This allows us to use from_atom37 to create the
        appropriate ProteinChain object with the appropriate atom37_mask.

        This function passes all kwargs to from_atom37.
        """
        if isinstance(backbone_atom_coordinates, torch.Tensor):
            backbone_atom_coordinates = backbone_atom_coordinates.cpu().numpy()
            if backbone_atom_coordinates.ndim == 4:
                if backbone_atom_coordinates.shape[0] != 1:
                    raise ValueError(
                        f"Cannot handle batched inputs, backbone_atom_coordinates has "
                        f"shape {backbone_atom_coordinates.shape}"
                    )
                backbone_atom_coordinates = backbone_atom_coordinates[0]

        assert isinstance(backbone_atom_coordinates, np.ndarray)
        assert backbone_atom_coordinates.ndim == 3
        assert backbone_atom_coordinates.shape[-2] == 3
        assert backbone_atom_coordinates.shape[-1] == 3

        atom37_positions = np.full(
            (backbone_atom_coordinates.shape[0], 37, 3),
            np.inf,
            dtype=backbone_atom_coordinates.dtype,
        )
        atom37_positions[:, :3, :] = backbone_atom_coordinates

        return cls.from_atom37(
            atom37_positions=atom37_positions,
            **kwargs,
        )

    @classmethod
    def from_pdb(
        cls,
        path: PathOrBuffer,
        chain_id: str = "detect",
        id: str | None = None,
        is_predicted: bool = False,
    ) -> "ProteinChain":
        """Return a ProteinChain object from an pdb file.

        Args:
            path (str | Path | io.TextIO): Path or buffer to read pdb file from. Should be uncompressed.
            id (str, optional): String identifier to assign to structure. Will attempt to infer otherwise.
            is_predicted (bool): If True, reads b factor as the confidence readout. Default: False.
            chain_id (str, optional): Select a chain corresponding to (author) chain id. "detect" uses the
                first detected chain
        """

        if id is not None:
            file_id = id
        else:
            match path:
                case Path() | str():
                    file_id = Path(path).with_suffix("").name
                case _:
                    file_id = "null"

        atom_array = PDBFile.read(path).get_structure(
            model=1, extra_fields=["b_factor"]
        )
        if chain_id == "detect":
            chain_id = atom_array.chain_id[0]
        atom_array = atom_array[
            bs.filter_amino_acids(atom_array)
            & ~atom_array.hetero
            & (atom_array.chain_id == chain_id)
        ]

        entity_id = 1  # Not supplied in PDBfiles

        sequence = "".join(
            (
                r
                if len(r := PDBData.protein_letters_3to1.get(monomer[0].res_name, "X"))
                == 1
                else "X"
            )
            for monomer in bs.residue_iter(atom_array)
        )
        num_res = len(sequence)

        atom_positions = np.full(
            [num_res, RC.atom_type_num, 3],
            np.nan,
            dtype=np.float32,
        )
        atom_mask = np.full(
            [num_res, RC.atom_type_num],
            False,
            dtype=bool,
        )
        residue_index = np.full([num_res], -1, dtype=np.int64)
        insertion_code = np.full([num_res], "", dtype="<U4")

        confidence = np.ones(
            [num_res],
            dtype=np.float32,
        )

        for i, res in enumerate(bs.residue_iter(atom_array)):
            chain = atom_array[atom_array.chain_id == chain_id]
            assert isinstance(chain, bs.AtomArray)

            res_index = res[0].res_id
            residue_index[i] = res_index
            insertion_code[i] = res[0].ins_code

            # Atom level features
            for atom in res:
                atom_name = atom.atom_name
                if atom_name == "SE" and atom.res_name == "MSE":
                    # Put the coords of the selenium atom in the sulphur column
                    atom_name = "SD"

                if atom_name in RC.atom_order:
                    atom_positions[i, RC.atom_order[atom_name]] = atom.coord
                    atom_mask[i, RC.atom_order[atom_name]] = True
                    if is_predicted and atom_name == "CA":
                        confidence[i] = atom.b_factor

        assert all(sequence), "Some residue name was not specified correctly"

        return cls(
            id=file_id,
            sequence=sequence,
            chain_id=chain_id,
            entity_id=entity_id,
            atom37_positions=atom_positions,
            atom37_mask=atom_mask,
            residue_index=residue_index,
            insertion_code=insertion_code,
            confidence=confidence,
        )

    @classmethod
    def from_rcsb(
        cls,
        pdb_id: str,
        chain_id: str = "detect",
    ):
        """Fetch a protein chain from the RCSB PDB database."""
        f: io.StringIO = rcsb.fetch(pdb_id, "pdb")  # type: ignore
        return cls.from_pdb(f, chain_id=chain_id, id=pdb_id)

    @classmethod
    def from_atomarray(
        cls,
        atom_array: bs.AtomArray,
        id: str | None = None,
    ) -> "ProteinChain":
        """A simple converter from bs.AtomArray -> ProteinChain.
        Uses PDB file format as intermediate."""
        atom_array = atom_array.copy()
        atom_array.box = None  # remove surrounding box, from_pdb won't handle this
        pdb_file = PDBFile()  # pyright: ignore
        pdb_file.set_structure(atom_array)

        buf = io.StringIO()
        pdb_file.write(buf)
        buf.seek(0)
        return cls.from_pdb(buf, id=id)

    def get_normalization_frame(self) -> Affine3D:
        """Given a set of coordinates, compute a single frame.
        Specifically, we compute the average position of the N, CA, and C atoms use those 3 points to construct a frame using the Gram-Schmidt algorithm. The average CA position is used as the origin of the frame.

        Returns:
            Affine3D: [] tensor of Affine3D frame
        """
        coords = torch.from_numpy(self.atom37_positions)
        frame = get_protein_normalization_frame(coords)

        return frame

    def apply_frame(self, frame: Affine3D) -> ProteinChain:
        """Given a frame, apply the frame to the protein's coordinates.

        Args:
            frame (Affine3D): [] tensor of Affine3D frame

        Returns:
            ProteinChain: Transformed protein chain
        """
        coords = torch.from_numpy(self.atom37_positions).to(frame.trans.dtype)
        coords = apply_frame_to_coords(coords, frame)
        atom37_positions = coords.numpy()
        return replace(self, atom37_positions=atom37_positions)

    def normalize_coordinates(self) -> ProteinChain:
        """Normalize the coordinates of the protein chain."""
        return self.apply_frame(self.get_normalization_frame())

    def infer_oxygen(self) -> ProteinChain:
        """Oxygen position is fixed given N, CA, C atoms. Infer it if not provided."""
        O_vector = torch.tensor([0.6240, -1.0613, 0.0103], dtype=torch.float32)
        N, CA, C = torch.from_numpy(self.atoms[["N", "CA", "C"]]).float().unbind(dim=1)
        N = torch.roll(N, -3)
        N[..., -1, :] = torch.nan

        # Get the frame defined by the CA-C-N atom
        frames = Affine3D.from_graham_schmidt(CA, C, N)
        O = frames.apply(O_vector)
        atom37_positions = self.atom37_positions.copy()
        atom37_mask = self.atom37_mask.copy()

        atom37_positions[:, RC.atom_order["O"]] = O.numpy()
        atom37_mask[:, RC.atom_order["O"]] = ~np.isnan(
            atom37_positions[:, RC.atom_order["O"]]
        ).any(-1)
        new_chain = replace(
            self, atom37_positions=atom37_positions, atom37_mask=atom37_mask
        )
        return new_chain

    @cached_property
    def inferred_cbeta(self) -> np.ndarray:
        """Infer cbeta positions based on N, C, CA."""
        N, CA, C = np.moveaxis(self.atoms[["N", "CA", "C"]], 1, 0)
        # See usage in trDesign codebase.
        # https://github.com/gjoni/trDesign/blob/f2d5930b472e77bfacc2f437b3966e7a708a8d37/02-GD/utils.py#L140
        CB = infer_CB(C, N, CA, 1.522, 1.927, -2.143)
        return CB

    def infer_cbeta(self, infer_cbeta_for_glycine: bool = False) -> ProteinChain:
        """Return a new chain with inferred CB atoms at all residues except GLY.

        Args:
            infer_cbeta_for_glycine (bool): If True, infers a beta carbon for glycine
                residues, even though that residue doesn't have one.  Default off.

                NOTE: The reason for having this switch in the first place
                is that sometimes we want a (inferred) CB coordinate for every residue,
                for example for making a pairwise distance matrix, or doing an RMSD
                calculation between two designs for a given structural template, w/
                CB atoms.
        """
        atom37_positions = self.atom37_positions.copy()
        atom37_mask = self.atom37_mask.copy()

        inferred_cbeta_positions = self.inferred_cbeta
        if not infer_cbeta_for_glycine:
            inferred_cbeta_positions[np.array(list(self.sequence)) == "G", :] = np.NAN

        atom37_positions[:, RC.atom_order["CB"]] = inferred_cbeta_positions
        atom37_mask[:, RC.atom_order["CB"]] = ~np.isnan(
            atom37_positions[:, RC.atom_order["CB"]]
        ).any(-1)
        new_chain = replace(
            self, atom37_positions=atom37_positions, atom37_mask=atom37_mask
        )
        return new_chain

    @cached_property
    def pdist_CA(self) -> np.ndarray:
        CA = self.atoms["CA"]
        pdist_CA = squareform(pdist(CA))
        return pdist_CA

    @cached_property
    def pdist_CB(self) -> np.ndarray:
        pdist_CB = squareform(pdist(self.inferred_cbeta))
        return pdist_CB

    @classmethod
    def as_complex(cls, chains: Sequence[ProteinChain]):
        raise RuntimeError(
            ".as_complex() has been deprecated in favor of .concat(). "
            ".concat() will eventually be deprecated in favor of ProteinComplex..."
        )

    @classmethod
    def concat(cls, chains: Sequence[ProteinChain]):
        def join_arrays(arrays: Sequence[np.ndarray], sep: np.ndarray):
            full_array = []
            for array in arrays:
                full_array.append(array)
                full_array.append(sep)
            full_array = full_array[:-1]
            return np.concatenate(full_array, 0)

        sep_tokens = {
            "residue_index": np.array([-1]),
            "insertion_code": np.array([""]),
            "atom37_positions": np.full([1, 37, 3], np.nan),
            "atom37_mask": np.zeros([1, 37]),
            "confidence": np.array([0]),
        }

        array_args: dict[str, np.ndarray] = {
            name: join_arrays([getattr(chain, name) for chain in chains], sep)
            for name, sep in sep_tokens.items()
        }

        return cls(
            id=chains[0].id,
            sequence=C.CHAIN_BREAK_STR.join(chain.sequence for chain in chains),
            chain_id="A",
            entity_id=None,
            **array_args,
        )

    def select_residue_indices(
        self, indices: list[int | str], ignore_x_mismatch: bool = False
    ) -> ProteinChain:
        numeric_indices = [
            idx if isinstance(idx, int) else int(idx[1:]) for idx in indices
        ]
        mask = np.isin(self.residue_index, numeric_indices)
        new = self[mask]
        mismatches = []
        for aa, idx in zip(new.sequence, indices):
            if isinstance(idx, int):
                continue
            if aa == "X" and ignore_x_mismatch:
                continue
            if aa != idx[0]:
                mismatches.append((aa, idx))
        if mismatches:
            mismatch_str = "; ".join(
                f"Position {idx[1:]}, Expected: {idx[0]}, Received: {aa}"
                for aa, idx in mismatches
            )
            raise RuntimeError(mismatch_str)

        return new
