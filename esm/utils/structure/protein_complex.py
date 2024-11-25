from __future__ import annotations

import io
import itertools
import re
import warnings
from dataclasses import asdict, dataclass, replace
from functools import cached_property
from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory
from typing import Any, Iterable, Sequence

import biotite.structure as bs
import brotli
import msgpack
import msgpack_numpy
import numpy as np
import torch
from biotite.database import rcsb
from biotite.structure.io.pdb import PDBFile

from esm.utils import residue_constants
from esm.utils.constants import esm3 as esm3_c
from esm.utils.misc import slice_python_object_as_numpy
from esm.utils.structure.affine3d import Affine3D
from esm.utils.structure.aligner import Aligner
from esm.utils.structure.metrics import (
    compute_gdt_ts,
    compute_lddt_ca,
)
from esm.utils.structure.protein_chain import (
    PathOrBuffer,
    ProteinChain,
)
from esm.utils.structure.protein_structure import (
    index_by_atom_name,
)

msgpack_numpy.patch()

SINGLE_LETTER_CHAIN_IDS = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
)


def protein_chain_to_protein_complex(chain: ProteinChain) -> ProteinComplex:
    if "|" not in chain.sequence:
        return ProteinComplex.from_chains([chain])
    chain_breaks = np.array(list(chain.sequence)) == "|"
    chain_break_inds = np.where(chain_breaks)[0]
    chain_break_inds = np.concatenate([[0], chain_break_inds, [len(chain)]])
    chain_break_inds = np.array(list(zip(chain_break_inds[:-1], chain_break_inds[1:])))
    complex_chains = []
    for start, end in chain_break_inds:
        if start != 0:
            start += 1
        complex_chains.append(chain[start:end])
    complex_chains = [
        ProteinChain.from_atom37(
            chain.atom37_positions,
            sequence=chain.sequence,
            chain_id=SINGLE_LETTER_CHAIN_IDS[i],
            entity_id=i,
        )
        for i, chain in enumerate(complex_chains)
    ]
    return ProteinComplex.from_chains(complex_chains)


@dataclass
class ProteinComplexMetadata:
    entity_lookup: dict[int, int]
    chain_lookup: dict[int, str]
    chain_boundaries: list[tuple[int, int]]


@dataclass
class DockQSingleScore:
    native_chains: tuple[str, str]
    DockQ: float
    interface_rms: float
    ligand_rms: float
    fnat: float
    fnonnat: float
    clashes: float
    F1: float
    DockQ_F1: float


@dataclass
class DockQResult:
    total_dockq: float
    native_interfaces: int
    chain_mapping: dict[str, str]
    interfaces: dict[tuple[str, str], DockQSingleScore]
    # zip(aligned.chain_iter(), native.chain_iter()) gives you the pairing
    # aligned.rmsd(native) should give you a low rmsd irrespective of shuffling
    aligned: ProteinComplex
    aligned_rmsd: float


class AtomIndexer:
    def __init__(self, structure: ProteinComplex, property: str, dim: int):
        self.structure = structure
        self.property = property
        self.dim = dim

    def __getitem__(self, atom_names: str | list[str]) -> np.ndarray:
        return index_by_atom_name(
            getattr(self.structure, self.property), atom_names, self.dim
        )


@dataclass
class ProteinComplex:
    """Dataclass with atom37 representation of an entire protein complex."""

    id: str
    sequence: str
    entity_id: np.ndarray  # entities map to unique sequences
    chain_id: np.ndarray  # multiple chains might share an entity id
    sym_id: np.ndarray  # complexes might be copies of the same chain
    residue_index: np.ndarray
    insertion_code: np.ndarray
    atom37_positions: np.ndarray
    atom37_mask: np.ndarray
    confidence: np.ndarray
    metadata: ProteinComplexMetadata

    def __post_init__(self):
        l = len(self.sequence)
        assert self.atom37_positions.shape[0] == l, (self.atom37_positions.shape, l)
        assert self.atom37_mask.shape[0] == l, (self.atom37_mask.shape, l)
        assert self.residue_index.shape[0] == l, (self.residue_index.shape, l)
        assert self.insertion_code.shape[0] == l, (self.insertion_code.shape, l)
        assert self.confidence.shape[0] == l, (self.confidence.shape, l)
        assert self.entity_id.shape[0] == l, (self.entity_id.shape, l)
        assert self.chain_id.shape[0] == l, (self.chain_id.shape, l)
        assert self.sym_id.shape[0] == l, (self.sym_id.shape, l)

    def __getitem__(self, idx: int | list[int] | slice | np.ndarray):
        """This function slices protein complexes without consideration of chain breaks
        NOTE: When slicing with a boolean mask, it's possible that the output array won't
        be the expected length. This is because we do our best to preserve chainbreak tokens.
        """

        if isinstance(idx, int):
            idx = [idx]
        if isinstance(idx, list):
            raise ValueError(
                "ProteinComplex doesn't supports indexing with lists of indices"
            )

        if isinstance(idx, np.ndarray):
            is_chainbreak = np.asarray([s == "|" for s in self.sequence])
            idx = idx.astype(bool) | is_chainbreak

        complex = self._unsafe_slice(idx)
        if len(complex) == 0:
            return complex

        # detect runs of chainbreaks by searching for instances of '||' in complex.sequence
        chainbreak_runs = np.asarray(
            [
                complex.sequence[i : i + 2] == "||"
                for i in range(len(complex.sequence) - 1)
            ]
            + [complex.sequence[-1] == "|"]
        )
        # We should remove as many chainbreaks as possible from the start of the sequence
        for i in range(len(chainbreak_runs)):
            if complex.sequence[i] == "|":
                chainbreak_runs[i] = True
            else:
                break
        complex = complex._unsafe_slice(~chainbreak_runs)
        return complex

    def _unsafe_slice(self, idx: int | list[int] | slice | np.ndarray):
        sequence = slice_python_object_as_numpy(self.sequence, idx)
        return replace(
            self,
            sequence=sequence,
            entity_id=self.entity_id[..., idx],
            chain_id=self.chain_id[..., idx],
            sym_id=self.sym_id[..., idx],
            residue_index=self.residue_index[..., idx],
            insertion_code=self.insertion_code[..., idx],
            atom37_positions=self.atom37_positions[..., idx, :, :],
            atom37_mask=self.atom37_mask[..., idx, :],
            confidence=self.confidence[..., idx],
        )

    def __len__(self):
        return len(self.sequence)

    @cached_property
    def atoms(self) -> AtomIndexer:
        return AtomIndexer(self, property="atom37_positions", dim=-2)

    def chain_iter(self) -> Iterable[ProteinChain]:
        boundaries = [i for i, s in enumerate(self.sequence) if s == "|"]
        boundaries = [-1, *boundaries, len(self)]
        for i in range(len(boundaries) - 1):
            c = self.__getitem__(slice(boundaries[i] + 1, boundaries[i + 1]))
            yield c.as_chain()

    def as_chain(self, force_conversion: bool = False) -> ProteinChain:
        """Convert the ProteinComplex to a ProteinChain.

        Args:
            force_conversion (bool): Forces the conversion into a protein chain even if the complex has multiple chains.
                The purpose of this is to use ProteinChain specific functions (like cbeta_contacts).

        """
        if not force_conversion:
            assert len(np.unique(self.chain_id)) == 1, f"{self.id}"
            assert len(np.unique(self.entity_id)) == 1, f"{self.id}"
            if self.chain_id[0] not in self.metadata.chain_lookup:
                warnings.warn("Chain ID not found in metadata, using 'A' as default")
            if self.entity_id[0] not in self.metadata.entity_lookup:
                warnings.warn("Entity ID not found in metadata, using None as default")
            chain_id = self.metadata.chain_lookup.get(self.chain_id[0], "A")
            entity_id = self.metadata.entity_lookup.get(self.entity_id[0], None)
        else:
            chain_id = "A"
            entity_id = None

        return ProteinChain(
            id=self.id,
            sequence=self.sequence,
            chain_id=chain_id,
            entity_id=entity_id,
            atom37_positions=self.atom37_positions,
            atom37_mask=self.atom37_mask,
            residue_index=self.residue_index,
            insertion_code=self.insertion_code,
            confidence=self.confidence,
        )

    @classmethod
    def from_pdb(cls, path: PathOrBuffer, id: str | None = None) -> "ProteinComplex":
        atom_array = PDBFile.read(path).get_structure(
            model=1, extra_fields=["b_factor"]
        )

        chains = []
        for chain in bs.chain_iter(atom_array):
            chain = chain[~chain.hetero]
            if len(chain) == 0:
                continue
            chains.append(ProteinChain.from_atomarray(chain, id))
        return ProteinComplex.from_chains(chains)

    @classmethod
    def from_rcsb(cls, pdb_id: str):
        """Fetch a protein complex from the RCSB PDB database."""
        f: io.StringIO = rcsb.fetch(pdb_id, "pdb")  # type: ignore
        return cls.from_pdb(f, id=pdb_id)

    def to_pdb(self, path: PathOrBuffer, include_insertions: bool = True):
        atom_array = None
        for chain in self.chain_iter():
            carr = (
                chain.atom_array
                if include_insertions
                else chain.atom_array_no_insertions
            )
            atom_array = carr if atom_array is None else atom_array + carr
        f = PDBFile()
        f.set_structure(atom_array)
        f.write(path)

    def to_pdb_string(self, include_insertions: bool = True) -> str:
        buf = io.StringIO()
        self.to_pdb(buf, include_insertions=include_insertions)
        buf.seek(0)
        return buf.read()

    def normalize_chain_ids_for_pdb(self):
        # Since PDB files have 1-letter chain IDs and don't support the idea of a symmetric index,
        # we can normalize it instead which might be necessary for DockQ and to_pdb.
        ids = SINGLE_LETTER_CHAIN_IDS
        chains = []
        for i, chain in enumerate(self.chain_iter()):
            chain.chain_id = ids[i]
            if i > len(ids):
                raise RuntimeError("Too many chains to write to PDB file")
            chains.append(chain)

        return ProteinComplex.from_chains(chains)

    def state_dict(self, backbone_only=False):
        """This state dict is optimized for storage, so it turns things to fp16 whenever
        possible. Note that we also only support int32 residue indices, I'm hoping we don't
        need more than 2**32 residues..."""
        dct = {k: v for k, v in vars(self).items()}
        for k, v in dct.items():
            if isinstance(v, np.ndarray):
                match v.dtype:
                    case np.int64:
                        dct[k] = v.astype(np.int32)
                    case np.float64 | np.float32:
                        dct[k] = v.astype(np.float16)
                    case _:
                        pass
            elif isinstance(v, ProteinComplexMetadata):
                dct[k] = asdict(v)
        dct["atom37_positions"] = dct["atom37_positions"][dct["atom37_mask"]]
        return dct

    def to_blob(self, backbone_only=False) -> bytes:
        return brotli.compress(msgpack.dumps(self.state_dict(backbone_only)), quality=5)

    @classmethod
    def from_state_dict(cls, dct):
        atom37 = np.full((*dct["atom37_mask"].shape, 3), np.nan)
        atom37[dct["atom37_mask"]] = dct["atom37_positions"]
        dct["atom37_positions"] = atom37
        dct = {
            k: (v.astype(np.float32) if k in ["atom37_positions", "confidence"] else v)
            for k, v in dct.items()
        }
        dct["metadata"] = ProteinComplexMetadata(**dct["metadata"])
        return cls(**dct)

    @classmethod
    def from_blob(cls, input: Path | str | io.BytesIO | bytes):
        """NOTE(@zlin): blob + sparse coding + brotli + fp16 reduces memory
        of chains from 52G/1M chains to 20G/1M chains, I think this is a good first
        shot at compressing and dumping chains to disk. I'm sure there's better ways."""
        match input:
            case Path() | str():
                bytes = Path(input).read_bytes()
            case io.BytesIO():
                bytes = input.getvalue()
            case _:
                bytes = input
        return cls.from_state_dict(
            msgpack.loads(brotli.decompress(bytes), strict_map_key=False)
        )

    @classmethod
    def from_chains(cls, chains: Sequence[ProteinChain]):
        if not chains:
            raise ValueError(
                "Cannot create a ProteinComplex from an empty list of chains"
            )

        # TODO: Make a proper protein complex class
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
            "atom37_mask": np.zeros([1, 37], dtype=bool),
            "confidence": np.array([0]),
        }

        array_args: dict[str, np.ndarray] = {
            name: join_arrays([getattr(chain, name) for chain in chains], sep)
            for name, sep in sep_tokens.items()
        }

        multimer_arrays = []
        chain2num_max = -1
        chain2num = {}
        ent2num_max = -1
        ent2num = {}
        total_index = 0
        chain_boundaries = []
        for i, c in enumerate(chains):
            num_res = c.residue_index.shape[0]
            if c.chain_id not in chain2num:
                chain2num[c.chain_id] = (chain2num_max := chain2num_max + 1)
            chain_id_array = np.full([num_res], chain2num[c.chain_id], dtype=np.int64)

            if c.entity_id is None:
                entity_num = (ent2num_max := ent2num_max + 1)
            else:
                if c.entity_id not in ent2num:
                    ent2num[c.entity_id] = (ent2num_max := ent2num_max + 1)
                entity_num = ent2num[c.entity_id]
            entity_id_array = np.full([num_res], entity_num, dtype=np.int64)

            sym_id_array = np.full([num_res], i, dtype=np.int64)

            multimer_arrays.append(
                {
                    "chain_id": chain_id_array,
                    "entity_id": entity_id_array,
                    "sym_id": sym_id_array,
                }
            )

            chain_boundaries.append((total_index, total_index + num_res))
            total_index += num_res + 1

        sep = np.array([-1])
        update = {
            name: join_arrays([dct[name] for dct in multimer_arrays], sep=sep)
            for name in ["chain_id", "entity_id", "sym_id"]
        }
        array_args.update(update)

        metadata = ProteinComplexMetadata(
            chain_boundaries=chain_boundaries,
            chain_lookup={v: k for k, v in chain2num.items()},
            entity_lookup={v: k for k, v in ent2num.items()},
        )

        return cls(
            id=chains[0].id,
            sequence=esm3_c.CHAIN_BREAK_STR.join(chain.sequence for chain in chains),
            metadata=metadata,
            **array_args,
        )

    def infer_oxygen(self) -> ProteinComplex:
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

        atom37_positions[:, residue_constants.atom_order["O"]] = O.numpy()
        atom37_mask[:, residue_constants.atom_order["O"]] = ~np.isnan(
            atom37_positions[:, residue_constants.atom_order["O"]]
        ).any(-1)
        new_chain = replace(
            self, atom37_positions=atom37_positions, atom37_mask=atom37_mask
        )
        return new_chain

    @classmethod
    def concat(cls, objs: list[ProteinComplex]) -> ProteinComplex:
        pdb_ids = [obj.id for obj in objs]
        if len(set(pdb_ids)) > 1:
            raise RuntimeError(
                "Concatention of protein complexes across different PDB ids is unsupported"
            )
        return ProteinComplex.from_chains(
            list(itertools.chain.from_iterable(obj.chain_iter() for obj in objs))
        )

    def _sanity_check_complexes_are_comparable(self, other: ProteinComplex):
        assert len(self) == len(other), "Protein complexes must have the same length"
        assert len(list(self.chain_iter())) == len(
            list(other.chain_iter())
        ), "Protein complexes must have the same number of chains"

    def lddt_ca(
        self,
        target: ProteinComplex,
        mobile_inds: list[int] | np.ndarray | None = None,
        target_inds: list[int] | np.ndarray | None = None,
        compute_chain_assignment: bool = True,
        **kwargs,
    ) -> float | np.ndarray:
        """Compute the LDDT between this protein complex and another.

        Arguments:
            target (ProteinComplex): The other protein complex to compare to.
            mobile_inds (list[int], np.ndarray, optional): The indices of the mobile atoms to align. These are NOT residue indices
            target_inds (list[int], np.ndarray, optional): The indices of the target atoms to align. These are NOT residue indices

        Returns:
            float | np.ndarray: The LDDT score between the two protein chains, either
                a single float or per-residue LDDT scores if `per_residue` is True.
        """
        if compute_chain_assignment:
            aligned = self.dockq(target).aligned
        else:
            aligned = self
        lddt = compute_lddt_ca(
            torch.tensor(aligned.atom37_positions[mobile_inds]).unsqueeze(0),
            torch.tensor(target.atom37_positions[target_inds]).unsqueeze(0),
            torch.tensor(aligned.atom37_mask[mobile_inds]).unsqueeze(0),
            **kwargs,
        )
        return float(lddt) if lddt.numel() == 1 else lddt.numpy().flatten()

    def gdt_ts(
        self,
        target: ProteinComplex,
        mobile_inds: list[int] | np.ndarray | None = None,
        target_inds: list[int] | np.ndarray | None = None,
        compute_chain_assignment: bool = True,
        **kwargs,
    ) -> float | np.ndarray:
        """Compute the GDT_TS between this protein complex and another.

        Arguments:
            target (ProteinComplex): The other protein complex to compare to.
            mobile_inds (list[int], np.ndarray, optional): The indices of the mobile atoms to align. These are NOT residue indices
            target_inds (list[int], np.ndarray, optional): The indices of the target atoms to align. These are NOT residue indices

        Returns:
            float: The GDT_TS score between the two protein chains.
        """
        if compute_chain_assignment:
            aligned = self.dockq(target).aligned
        else:
            aligned = self
        gdt_ts = compute_gdt_ts(
            mobile=torch.tensor(
                index_by_atom_name(aligned.atom37_positions[mobile_inds], "CA"),
                dtype=torch.float32,
            ).unsqueeze(0),
            target=torch.tensor(
                index_by_atom_name(target.atom37_positions[target_inds], "CA"),
                dtype=torch.float32,
            ).unsqueeze(0),
            atom_exists_mask=torch.tensor(
                index_by_atom_name(aligned.atom37_mask[mobile_inds], "CA", dim=-1)
                & index_by_atom_name(target.atom37_mask[target_inds], "CA", dim=-1)
            ).unsqueeze(0),
            **kwargs,
        )
        return float(gdt_ts) if gdt_ts.numel() == 1 else gdt_ts.numpy().flatten()

    def dockq(self, native: ProteinComplex):
        # This function uses dockqv2 to compute the DockQ score. Because it does a mapping
        # over all possible chains, it's quite slow. Be careful not to use this in an inference loop
        # or something that requires fast scoring. It defaults to 8 CPUs.

        try:
            pass
        except BaseException:
            raise RuntimeError(
                "DockQ is not installed. Please update your environment."
            )
        self._sanity_check_complexes_are_comparable(native)

        def sanity_check_chain_ids(pc: ProteinComplex):
            ids = []
            for i, chain in enumerate(pc.chain_iter()):
                if i > len(SINGLE_LETTER_CHAIN_IDS):
                    raise ValueError("Too many chains to write to PDB file")
                if len(chain.chain_id) > 1:
                    raise ValueError(
                        "We only supports single letter chain IDs for DockQ"
                    )
                ids.append(chain.chain_id)
            if len(set(ids)) != len(ids):
                raise ValueError(f"Duplicate chain IDs in protein complex: {ids}")
            return ids

        sanity_check_chain_ids(self)
        sanity_check_chain_ids(native)

        with TemporaryDirectory() as tdir:
            dir = Path(tdir)
            self.to_pdb(dir / "self.pdb")
            native.to_pdb(dir / "native.pdb")

            output = check_output(["DockQ", dir / "self.pdb", dir / "native.pdb"])
        lines = output.decode().split("\n")

        # Remove the header comments
        start_index = next(
            i for i, line in enumerate(lines) if line.startswith("Model")
        )
        lines = lines[start_index:]

        result = {}
        interfaces = []
        current_interface: dict = {}

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Model  :"):
                pass  # Tmp pdb file location, it's useless...
            elif line.startswith("Native :"):
                pass  # Tmp pdb file location, it's useless...
            elif line.startswith("Total DockQ"):
                total_dockq_match = re.search(
                    r"Total DockQ over (\d+) native interfaces: ([\d.]+) with (.*) model:native mapping",
                    line,
                )
                if total_dockq_match:
                    result["value"] = float(total_dockq_match.group(2))
                    result["native interfaces"] = int(total_dockq_match.group(1))
                    native_chains, self_chains = total_dockq_match.group(3).split(":")
                    result["mapping"] = dict(zip(native_chains, self_chains))
                else:
                    raise RuntimeError(
                        "Failed to parse DockQ output, maybe your DockQ version is wrong?"
                    )
            elif line.startswith("Native chains:"):
                if current_interface:
                    interfaces.append(current_interface)
                current_interface = {
                    "Native chains": line.split(":")[1].strip().split(", ")
                }
            elif line.startswith("Model chains:"):
                current_interface["Model chains"] = (
                    line.split(":")[1].strip().split(", ")
                )
            elif ":" in line:
                key, value = line.split(":", 1)
                current_interface[key.strip()] = float(value.strip())

        if current_interface:
            interfaces.append(current_interface)

        def parse_dict(d: dict[str, Any]) -> DockQSingleScore:
            return DockQSingleScore(
                native_chains=tuple(d["Native chains"]),  # type: ignore
                DockQ=float(d["DockQ"]),
                interface_rms=float(d["irms"]),
                ligand_rms=float(d["Lrms"]),  # Note the capitalization difference
                fnat=float(d["fnat"]),
                fnonnat=float(d["fnonnat"]),
                clashes=float(d["clashes"]),
                F1=float(d["F1"]),
                DockQ_F1=float(d["DockQ_F1"]),
            )

        inv_mapping = {v: k for k, v in result["mapping"].items()}

        self_chain_map = {c.chain_id: c for c in self.chain_iter()}
        realigned = []
        for chain in native.chain_iter():
            realigned.append(self_chain_map[inv_mapping[chain.chain_id]])

        realigned = ProteinComplex.from_chains(realigned)
        aligner = Aligner(realigned, native)
        realigned = aligner.apply(realigned)

        result = DockQResult(
            total_dockq=result["value"],
            native_interfaces=result["native interfaces"],
            chain_mapping=result["mapping"],
            interfaces={
                (i["Model chains"][0], i["Model chains"][1]): parse_dict(i)
                for i in interfaces
            },
            aligned=realigned,
            aligned_rmsd=aligner.rmsd,
        )

        return result
