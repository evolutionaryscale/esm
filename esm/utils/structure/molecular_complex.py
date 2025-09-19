from __future__ import annotations

import io
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any, List

import biotite.structure.io.pdbx as pdbx
import brotli
import msgpack
import numpy as np
import torch

from esm.utils import residue_constants
from esm.utils.structure.metrics import compute_lddt, compute_rmsd
from esm.utils.structure.protein_complex import ProteinComplex, ProteinComplexMetadata


@dataclass
class MolecularComplexResult:
    """Result of molecular complex folding"""

    complex: MolecularComplex
    plddt: torch.Tensor | None = None
    ptm: float | None = None
    iptm: float | None = None
    pae: torch.Tensor | None = None
    distogram: torch.Tensor | None = None
    pair_chains_iptm: torch.Tensor | None = None
    output_embedding_sequence: torch.Tensor | None = None
    output_embedding_pair_pooled: torch.Tensor | None = None


@dataclass
class MolecularComplexMetadata:
    """Metadata for MolecularComplex objects."""

    entity_lookup: dict[int, str]
    chain_lookup: dict[int, str]
    assembly_composition: dict[str, list[str]] | None = None


@dataclass
class Molecule:
    """Represents a single molecule/token within a MolecularComplex."""

    token: str
    token_idx: int
    atom_positions: np.ndarray  # [N_atoms, 3]
    atom_elements: np.ndarray  # [N_atoms] element strings
    residue_type: int
    molecule_type: int  # PROTEIN=0, RNA=1, DNA=2, LIGAND=3
    confidence: float


@dataclass(frozen=True)
class MolecularComplex:
    """
    Dataclass representing a molecular complex with support for proteins, nucleic acids, and ligands.

    Uses a flat atom representation with token-based sequence indexing, supporting all atom types
    beyond the traditional atom37 protein representation.
    """

    id: str
    sequence: List[str]  # Token sequence like ['MET', 'LYS', 'A', 'G', 'ATP']

    # Flat atom arrays - simplified representation
    atom_positions: np.ndarray  # [N_atoms, 3] 3D coordinates
    atom_elements: np.ndarray  # [N_atoms] element strings

    # Token-to-atom mapping for efficient access
    token_to_atoms: np.ndarray  # [N_tokens, 2] start/end indices into atoms array

    # Confidence data
    plddt: np.ndarray  # Per-token confidence scores [N_tokens]

    # Metadata
    metadata: MolecularComplexMetadata

    def __post_init__(self):
        """Validate array dimensions."""
        n_tokens = len(self.sequence)
        assert (
            self.token_to_atoms.shape[0] == n_tokens
        ), f"token_to_atoms shape {self.token_to_atoms.shape} != {n_tokens} tokens"
        assert (
            self.plddt.shape[0] == n_tokens
        ), f"plddt shape {self.plddt.shape} != {n_tokens} tokens"

    def __len__(self) -> int:
        """Return number of tokens."""
        return len(self.sequence)

    def __getitem__(self, idx: int) -> Molecule:
        """Access individual molecules/tokens by index."""
        if idx >= len(self.sequence) or idx < 0:
            raise IndexError(
                f"Token index {idx} out of range for {len(self.sequence)} tokens"
            )

        token = self.sequence[idx]
        start_atom, end_atom = self.token_to_atoms[idx]

        # Extract atom data for this token
        token_atom_positions = self.atom_positions[start_atom:end_atom]
        token_atom_elements = self.atom_elements[start_atom:end_atom]

        # Default values for residue/molecule type (would be extended based on actual implementation)
        residue_type = 0  # Default to standard residue
        molecule_type = 0  # Default to protein

        return Molecule(
            token=token,
            token_idx=idx,
            atom_positions=token_atom_positions,
            atom_elements=token_atom_elements,
            residue_type=residue_type,
            molecule_type=molecule_type,
            confidence=self.plddt[idx],
        )

    @property
    def atom_coordinates(self) -> np.ndarray:
        """Get flat array of all atom coordinates [N_atoms, 3]."""
        return self.atom_positions

    # Conversion methods
    @classmethod
    def from_protein_complex(cls, pc: ProteinComplex) -> "MolecularComplex":
        """Convert a ProteinComplex to MolecularComplex.

        Args:
            pc: ProteinComplex object with atom37 representation

        Returns:
            MolecularComplex with flat atom arrays and token-based indexing
        """
        from esm.utils import residue_constants

        # Extract sequence without chain breaks
        sequence_no_breaks = pc.sequence.replace("|", "")
        sequence_tokens = [
            residue_constants.restype_1to3.get(aa, "UNK") for aa in sequence_no_breaks
        ]

        # Convert atom37 to flat arrays
        flat_positions = []
        flat_elements = []
        token_to_atoms = []

        atom_idx = 0
        residue_idx = 0

        for i, aa in enumerate(pc.sequence):
            if aa == "|":
                # Skip chain break tokens
                continue

            # Get atom37 positions and mask for this residue
            res_positions = pc.atom37_positions[residue_idx]  # [37, 3]
            res_mask = pc.atom37_mask[residue_idx]  # [37]

            # Track start position for this token
            token_start = atom_idx

            # Process each atom type in atom37 representation
            for atom_type_idx, atom_name in enumerate(residue_constants.atom_types):
                if res_mask[atom_type_idx]:  # Atom is present
                    # Add position
                    flat_positions.append(res_positions[atom_type_idx])

                    # Determine element from atom name
                    element = (
                        atom_name[0] if atom_name else "C"
                    )  # First character is element
                    flat_elements.append(element)

                    atom_idx += 1

            # Record token-to-atom mapping [start_idx, end_idx)
            token_to_atoms.append([token_start, atom_idx])
            residue_idx += 1

        # Convert to numpy arrays
        atom_positions = np.array(flat_positions, dtype=np.float32)
        atom_elements = np.array(flat_elements, dtype=object)
        token_to_atoms_array = np.array(token_to_atoms, dtype=np.int32)

        # Extract confidence scores (skip chain breaks)
        confidence_scores = []
        residue_idx = 0
        for aa in pc.sequence:
            if aa != "|":
                confidence_scores.append(pc.confidence[residue_idx])
                residue_idx += 1

        confidence_array = np.array(confidence_scores, dtype=np.float32)

        # Create metadata - convert entity IDs to strings for MolecularComplexMetadata
        entity_lookup_str = {k: str(v) for k, v in pc.metadata.entity_lookup.items()}
        metadata = MolecularComplexMetadata(
            entity_lookup=entity_lookup_str,
            chain_lookup=pc.metadata.chain_lookup,
            assembly_composition=pc.metadata.assembly_composition,
        )

        return cls(
            id=pc.id,
            sequence=sequence_tokens,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            token_to_atoms=token_to_atoms_array,
            plddt=confidence_array,
            metadata=metadata,
        )

    def to_protein_complex(self) -> ProteinComplex:
        """Convert MolecularComplex back to ProteinComplex format.

        Extracts only protein tokens and converts from flat atom representation
        back to atom37 format used by ProteinComplex.

        Returns:
            ProteinComplex with protein residues only, excluding ligands/nucleic acids
        """
        from esm.utils import residue_constants

        # No need for element mapping - already using element characters

        # Filter for protein tokens only (skip ligands, nucleic acids)
        protein_tokens = []
        protein_indices = []

        for i, token in enumerate(self.sequence):
            # Check if token is a standard 3-letter amino acid code
            if token in residue_constants.restype_3to1:
                protein_tokens.append(token)
                protein_indices.append(i)

        if not protein_tokens:
            raise ValueError("No protein tokens found in MolecularComplex")

        n_residues = len(protein_tokens)

        # Initialize atom37 arrays
        atom37_positions = np.full((n_residues, 37, 3), np.nan, dtype=np.float32)
        atom37_mask = np.zeros((n_residues, 37), dtype=bool)

        # Convert tokens back to single-letter sequence
        single_letter_sequence = "".join(
            [residue_constants.restype_3to1[token] for token in protein_tokens]
        )

        # Extract confidence scores for protein residues only
        protein_confidence = self.plddt[protein_indices]

        # Convert flat atoms back to atom37 representation
        for res_idx, token_idx in enumerate(protein_indices):
            token = self.sequence[token_idx]
            start_atom, end_atom = self.token_to_atoms[token_idx]

            # Get atom data for this residue
            res_atom_positions = self.atom_positions[start_atom:end_atom]

            # Reconstruct atom37 representation by exactly reversing the forward conversion logic
            # In from_protein_complex, atoms are added in atom_types order if present in mask
            # So we need to reconstruct the mask and positions in the same order
            atom_count = 0
            for atom_type_idx, atom_name in enumerate(residue_constants.atom_types):
                # Check if this atom type exists for this residue and was present
                residue_atoms = residue_constants.residue_atoms.get(token, [])
                if atom_name in residue_atoms:
                    # This atom type exists for this residue, so it should have been included
                    if atom_count < len(res_atom_positions):
                        atom37_positions[res_idx, atom_type_idx] = res_atom_positions[
                            atom_count
                        ]
                        atom37_mask[res_idx, atom_type_idx] = True
                        atom_count += 1

        # Create other required arrays for ProteinComplex
        # For simplicity, assume all protein residues belong to the same entity/chain
        entity_id = np.zeros(n_residues, dtype=np.int64)
        chain_id = np.zeros(n_residues, dtype=np.int64)
        sym_id = np.zeros(n_residues, dtype=np.int64)
        residue_index = np.arange(1, n_residues + 1, dtype=np.int64)
        insertion_code = np.array([""] * n_residues, dtype=object)

        # Create simplified protein complex metadata
        # Map the first entity/chain from molecular complex metadata
        protein_metadata = ProteinComplexMetadata(
            entity_lookup={0: 1},  # Single entity (int for ProteinComplexMetadata)
            chain_lookup={0: "A"},  # Single chain
            assembly_composition=self.metadata.assembly_composition,
        )

        return ProteinComplex(
            id=self.id,
            sequence=single_letter_sequence,
            entity_id=entity_id,
            chain_id=chain_id,
            sym_id=sym_id,
            residue_index=residue_index,
            insertion_code=insertion_code,
            atom37_positions=atom37_positions,
            atom37_mask=atom37_mask,
            confidence=protein_confidence,
            metadata=protein_metadata,
        )

    @classmethod
    def from_mmcif(cls, inp: str, id: str | None = None) -> "MolecularComplex":
        """Read MolecularComplex from mmcif file or string.

        Args:
            inp: Path to mmCIF file or mmCIF content as string
            id: Optional identifier to assign to the complex

        Returns:
            MolecularComplex with all molecules (proteins, ligands, nucleic acids)
        """
        from io import StringIO

        # Check if input is a file path or mmCIF string content
        if os.path.exists(inp):
            # Input is a file path
            mmcif_file = pdbx.CIFFile.read(inp)
        else:
            # Input is mmCIF string content
            mmcif_file = pdbx.CIFFile.read(StringIO(inp))

        # Get structure - handle missing model information gracefully
        try:
            structure = pdbx.get_structure(mmcif_file, model=1)
        except (KeyError, ValueError):
            # Fallback for mmCIF files without model information
            try:
                structure = pdbx.get_structure(mmcif_file)
            except Exception:
                # Last resort: use the first available model or all atoms
                structure = pdbx.get_structure(mmcif_file, model=None)
        # Type hint for pyright - structure is an AtomArray which is iterable
        if TYPE_CHECKING:
            structure: Any = structure

        # Get entity information from mmCIF
        entity_info = {}
        try:
            # Access the first block in CIFFile
            block = mmcif_file[0]
            if "entity" in block:
                entity_category = block["entity"]
                if "id" in entity_category and "type" in entity_category:
                    entity_ids = entity_category["id"]
                    entity_types = entity_category["type"]
                    # Convert CIFColumn to list for iteration
                    if hasattr(entity_ids, "__iter__") and hasattr(
                        entity_types, "__iter__"
                    ):
                        # Type annotation to help pyright understand these are iterable
                        entity_ids_list = list(entity_ids)  # type: ignore
                        entity_types_list = list(entity_types)  # type: ignore
                        for eid, etype in zip(entity_ids_list, entity_types_list):
                            entity_info[eid] = etype
        except Exception:
            pass

        # Initialize arrays for flat atom representation
        sequence_tokens = []
        flat_positions = []
        flat_elements = []
        token_to_atoms = []
        confidence_scores = []

        atom_idx = 0

        # Group atoms by chain and residue
        chain_residue_groups = {}
        for atom in structure:
            chain_id = atom.chain_id
            res_id = atom.res_id
            res_name = atom.res_name

            if chain_id not in chain_residue_groups:
                chain_residue_groups[chain_id] = {}
            if res_id not in chain_residue_groups[chain_id]:
                chain_residue_groups[chain_id][res_id] = {
                    "atoms": [],
                    "res_name": res_name,
                    "is_hetero": atom.hetero,
                }
            chain_residue_groups[chain_id][res_id]["atoms"].append(atom)

        # Process each chain and residue
        for chain_id in sorted(chain_residue_groups.keys()):
            residues = chain_residue_groups[chain_id]

            for res_id in sorted(residues.keys()):
                residue_data = residues[res_id]
                res_name = residue_data["res_name"]
                atoms = residue_data["atoms"]
                is_hetero = residue_data["is_hetero"]

                # Skip water molecules
                if res_name == "HOH":
                    continue

                # Determine token name
                if not is_hetero and res_name in residue_constants.restype_3to1:
                    # Standard amino acid
                    token_name = res_name
                elif res_name in ["A", "T", "G", "C", "U", "DA", "DT", "DG", "DC"]:
                    # Nucleotide
                    token_name = res_name
                else:
                    # Ligand or other molecule
                    token_name = res_name

                sequence_tokens.append(token_name)
                token_start = atom_idx

                # Add all atoms from this residue
                for atom in atoms:
                    flat_positions.append(atom.coord)

                    # Get element character
                    element = atom.element
                    flat_elements.append(element)

                    atom_idx += 1

                # Record token-to-atom mapping
                token_to_atoms.append([token_start, atom_idx])

                # Add confidence score (B-factor if available, otherwise 1.0)
                bfactor = getattr(atoms[0], "b_factor", 50.0) if atoms else 50.0
                confidence_scores.append(min(bfactor / 100.0, 1.0))

        # Convert to numpy arrays
        if not flat_positions:
            # Create minimal arrays if no atoms found
            atom_positions = np.zeros((0, 3), dtype=np.float32)
            atom_elements = np.zeros(0, dtype=object)
            token_to_atoms_array = np.zeros((len(sequence_tokens), 2), dtype=np.int32)
        else:
            atom_positions = np.array(flat_positions, dtype=np.float32)
            atom_elements = np.array(flat_elements, dtype=object)
            token_to_atoms_array = np.array(token_to_atoms, dtype=np.int32)

        confidence_array = np.array(confidence_scores, dtype=np.float32)

        # Create metadata
        metadata = MolecularComplexMetadata(
            entity_lookup=entity_info,
            chain_lookup={
                i: chain_id for i, chain_id in enumerate(chain_residue_groups.keys())
            },
            assembly_composition=None,
        )

        # Set complex ID - if input was a path, use the stem; otherwise use default
        if os.path.exists(inp):
            complex_id = id or Path(inp).stem
        else:
            complex_id = id or "complex_from_string"

        return cls(
            id=complex_id,
            sequence=sequence_tokens,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            token_to_atoms=token_to_atoms_array,
            plddt=confidence_array,
            metadata=metadata,
        )

    def to_mmcif(self) -> str:
        """Write MolecularComplex to mmcif string.

        Returns:
            String representation of the complex in mmCIF format
        """
        # No need for element mapping - already using element characters

        lines = []

        # Header
        lines.append(f"data_{self.id}")
        lines.append("#")
        lines.append(f"_entry.id    {self.id}")
        lines.append("#")

        # Structure metadata
        lines.append("_struct.entry_id                    {}".format(self.id))
        lines.append("_struct.title                       'Protein Structure'")
        lines.append("#")

        # Entity information
        entity_id = 1
        chain_counter = 0
        lines.append("loop_")
        lines.append("_entity.id")
        lines.append("_entity.type")
        lines.append("_entity.pdbx_description")

        # Determine entities based on sequence
        protein_tokens = []
        other_tokens = []

        for i, token in enumerate(self.sequence):
            if token in residue_constants.restype_3to1:
                protein_tokens.append((i, token))
            else:
                other_tokens.append((i, token))

        if protein_tokens:
            lines.append(f"{entity_id}   polymer   'Protein chain'")
            entity_id += 1

        for token in set(token for _, token in other_tokens):
            lines.append(f"{entity_id}   non-polymer   'Ligand {token}'")
            entity_id += 1

        lines.append("#")

        # Chain assignments
        lines.append("loop_")
        lines.append("_struct_asym.id")
        lines.append("_struct_asym.entity_id")

        chain_id = "A"
        if protein_tokens:
            lines.append(f"{chain_id}   1")
            chain_counter += 1
            chain_id = chr(ord(chain_id) + 1)

        entity_id = 2
        for token in set(token for _, token in other_tokens):
            lines.append(f"{chain_id}   {entity_id}")
            entity_id += 1
            chain_counter += 1
            if chain_counter < 26:
                chain_id = chr(ord(chain_id) + 1)

        lines.append("#")

        # Atom site information
        lines.append("loop_")
        lines.append("_atom_site.group_PDB")
        lines.append("_atom_site.id")
        lines.append("_atom_site.type_symbol")
        lines.append("_atom_site.label_atom_id")
        lines.append("_atom_site.label_alt_id")
        lines.append("_atom_site.label_comp_id")
        lines.append("_atom_site.label_asym_id")
        lines.append("_atom_site.label_entity_id")
        lines.append("_atom_site.label_seq_id")
        lines.append("_atom_site.pdbx_PDB_ins_code")
        lines.append("_atom_site.Cartn_x")
        lines.append("_atom_site.Cartn_y")
        lines.append("_atom_site.Cartn_z")
        lines.append("_atom_site.occupancy")
        lines.append("_atom_site.B_iso_or_equiv")
        lines.append("_atom_site.pdbx_PDB_model_num")
        lines.append("_atom_site.auth_seq_id")
        lines.append("_atom_site.auth_comp_id")
        lines.append("_atom_site.auth_asym_id")
        lines.append("_atom_site.auth_atom_id")

        atom_id = 1
        seq_id = 1
        chain_id = "A"
        entity_id = 1

        for token_idx, token in enumerate(self.sequence):
            start_atom, end_atom = self.token_to_atoms[token_idx]

            # Determine if this is a protein residue or ligand
            is_protein = token in residue_constants.restype_3to1
            group_pdb = "ATOM" if is_protein else "HETATM"
            current_entity_id = 1 if is_protein else 2  # Simplified entity assignment
            current_chain_id = "A" if is_protein else "B"  # Simplified chain assignment

            # Create atom names for this token
            atom_names = []
            if is_protein:
                # Use standard protein atom names
                res_atoms = residue_constants.residue_atoms.get(
                    token, ["N", "CA", "C", "O"]
                )
                atom_names = res_atoms[: end_atom - start_atom]
            else:
                # Generate generic atom names for ligands
                for i in range(end_atom - start_atom):
                    atom_names.append(f"C{i+1}")

            # Pad atom names if needed
            while len(atom_names) < (end_atom - start_atom):
                atom_names.append(f"X{len(atom_names)+1}")

            # Write atoms for this token
            for atom_idx_in_token, global_atom_idx in enumerate(
                range(start_atom, end_atom)
            ):
                pos = self.atom_positions[global_atom_idx]
                element_char = self.atom_elements[global_atom_idx]
                element_symbol = element_char if isinstance(element_char, str) else "C"

                atom_name = (
                    atom_names[atom_idx_in_token]
                    if atom_idx_in_token < len(atom_names)
                    else f"X{atom_idx_in_token+1}"
                )

                # Format atom site line
                bfactor = (
                    self.plddt[token_idx] * 100.0
                    if len(self.plddt) > token_idx
                    else 50.0
                )

                line = (
                    f"{group_pdb:<6} {atom_id:>5} {element_symbol:<2} {atom_name:<4} . "
                    f"{token:<3} {current_chain_id} {current_entity_id} {seq_id:>3} ? "
                    f"{pos[0]:>8.3f} {pos[1]:>8.3f} {pos[2]:>8.3f} 1.00 {bfactor:>6.2f} 1 "
                    f"{seq_id:>3} {token:<3} {current_chain_id} {atom_name:<4}"
                )
                lines.append(line)
                atom_id += 1

            seq_id += 1

        lines.append("#")
        return "\n".join(lines)

    def dockq(self, native: "MolecularComplex") -> Any:
        """Compute DockQ score against native structure.

        Args:
            native: Native MolecularComplex to compute DockQ against

        Returns:
            DockQ result containing score and alignment information
        """
        # Imports moved to top of file

        # Convert both complexes to ProteinComplex format for DockQ computation
        # This extracts only the protein portion and converts to PDB format
        try:
            self_pc = self.to_protein_complex()
            native_pc = native.to_protein_complex()
        except ValueError as e:
            raise ValueError(
                f"Cannot convert MolecularComplex to ProteinComplex for DockQ: {e}"
            )

        # Normalize chain IDs for PDB compatibility
        self_pc = self_pc.normalize_chain_ids_for_pdb()
        native_pc = native_pc.normalize_chain_ids_for_pdb()

        # Use the existing ProteinComplex.dockq() method
        try:
            dockq_result = self_pc.dockq(native_pc)
            return dockq_result
        except Exception:
            # Fallback to manual DockQ computation if ProteinComplex.dockq() fails
            return self._compute_dockq_manual(native)

    def _compute_dockq_manual(self, native: "MolecularComplex") -> Any:
        """Manual DockQ computation fallback."""
        # Imports moved to top of file

        # Convert both complexes to ProteinComplex format
        try:
            self_pc = self.to_protein_complex()
            native_pc = native.to_protein_complex()
        except ValueError as e:
            raise ValueError(
                f"Cannot convert MolecularComplex to ProteinComplex for DockQ: {e}"
            )

        # Normalize chain IDs for PDB compatibility
        self_pc = self_pc.normalize_chain_ids_for_pdb()
        native_pc = native_pc.normalize_chain_ids_for_pdb()

        # Write temporary PDB files and run DockQ
        with TemporaryDirectory() as tdir:
            dir_path = Path(tdir)
            self_pdb = dir_path / "self.pdb"
            native_pdb = dir_path / "native.pdb"

            # Write PDB files
            self_pc.to_pdb(self_pdb)
            native_pc.to_pdb(native_pdb)

            # Run DockQ
            try:
                output = check_output(["DockQ", str(self_pdb), str(native_pdb)])
                output_text = output.decode()

                # Parse DockQ output
                lines = output_text.split("\n")

                # Find the total DockQ score
                dockq_score = None
                for line in lines:
                    if "Total DockQ" in line:
                        match = re.search(r"Total DockQ.*: ([\d.]+)", line)
                        if match:
                            dockq_score = float(match.group(1))
                            break

                if dockq_score is None:
                    # Try to find individual DockQ scores
                    for line in lines:
                        if line.startswith("DockQ") and ":" in line:
                            try:
                                dockq_score = float(line.split(":")[1].strip())
                                break
                            except (ValueError, IndexError):
                                continue

                if dockq_score is None:
                    raise ValueError("Could not parse DockQ score from output")

                # Return a simple result structure
                return {
                    "total_dockq": dockq_score,
                    "raw_output": output_text,
                    "aligned": self,  # Return self as aligned structure
                }

            except FileNotFoundError:
                raise RuntimeError(
                    "DockQ is not installed. Please install DockQ to use this method."
                )
            except Exception as e:
                raise RuntimeError(f"DockQ computation failed: {e}")

    def rmsd(self, target: "MolecularComplex", **kwargs) -> float:
        """Compute RMSD against target structure.

        Args:
            target: Target MolecularComplex to compute RMSD against
            **kwargs: Additional arguments passed to compute_rmsd

        Returns:
            float: RMSD value between the two structures
        """
        # Imports moved to top of file

        # Ensure both complexes have the same number of tokens
        if len(self) != len(target):
            raise ValueError(
                f"Complexes must have the same number of tokens: {len(self)} vs {len(target)}"
            )

        # Extract center positions for each token (using centroid of atoms)
        mobile_coords = []
        target_coords = []
        atom_mask = []

        for i in range(len(self)):
            # Get atom positions for this token
            mobile_start, mobile_end = self.token_to_atoms[i]
            target_start, target_end = target.token_to_atoms[i]

            # Extract atom positions
            mobile_atoms = self.atom_positions[mobile_start:mobile_end]
            target_atoms = target.atom_positions[target_start:target_end]

            # Check if both tokens have atoms
            if len(mobile_atoms) == 0 or len(target_atoms) == 0:
                # Skip tokens with no atoms
                continue

            # For simplicity, use the centroid of atoms as the representative position
            mobile_center = mobile_atoms.mean(axis=0)
            target_center = target_atoms.mean(axis=0)

            mobile_coords.append(mobile_center)
            target_coords.append(target_center)
            atom_mask.append(True)

        if len(mobile_coords) == 0:
            raise ValueError("No valid atoms found for RMSD computation")

        # Convert to tensors
        mobile_tensor = torch.from_numpy(np.stack(mobile_coords, axis=0)).unsqueeze(
            0
        )  # [1, N, 3]
        target_tensor = torch.from_numpy(np.stack(target_coords, axis=0)).unsqueeze(
            0
        )  # [1, N, 3]
        mask_tensor = torch.tensor(atom_mask, dtype=torch.bool).unsqueeze(0)  # [1, N]

        # Compute RMSD using existing infrastructure
        rmsd_value = compute_rmsd(
            mobile=mobile_tensor,
            target=target_tensor,
            atom_exists_mask=mask_tensor,
            reduction="batch",
            **kwargs,
        )

        return float(rmsd_value)

    def lddt_ca(self, target: "MolecularComplex", **kwargs) -> float:
        """Compute LDDT score against target structure.

        Args:
            target: Target MolecularComplex to compute LDDT against
            **kwargs: Additional arguments passed to compute_lddt

        Returns:
            float: LDDT value between the two structures
        """
        # Imports moved to top of file

        # Ensure both complexes have the same number of tokens
        if len(self) != len(target):
            raise ValueError(
                f"Complexes must have the same number of tokens: {len(self)} vs {len(target)}"
            )

        # Extract center positions for each token (using centroid of atoms)
        mobile_coords = []
        target_coords = []
        atom_mask = []

        for i in range(len(self)):
            # Get atom positions for this token
            mobile_start, mobile_end = self.token_to_atoms[i]
            target_start, target_end = target.token_to_atoms[i]

            # Extract atom positions
            mobile_atoms = self.atom_positions[mobile_start:mobile_end]
            target_atoms = target.atom_positions[target_start:target_end]

            # Check if both tokens have atoms
            if len(mobile_atoms) == 0 or len(target_atoms) == 0:
                # Skip tokens with no atoms
                mobile_coords.append(np.full(3, np.nan))
                target_coords.append(np.full(3, np.nan))
                atom_mask.append(False)
                continue

            # For simplicity, use the centroid of atoms as the representative position
            mobile_center = mobile_atoms.mean(axis=0)
            target_center = target_atoms.mean(axis=0)

            mobile_coords.append(mobile_center)
            target_coords.append(target_center)
            atom_mask.append(True)

        if not any(atom_mask):
            raise ValueError("No valid atoms found for LDDT computation")

        # Convert to tensors
        mobile_tensor = torch.from_numpy(np.stack(mobile_coords, axis=0)).unsqueeze(
            0
        )  # [1, N, 3]
        target_tensor = torch.from_numpy(np.stack(target_coords, axis=0)).unsqueeze(
            0
        )  # [1, N, 3]
        mask_tensor = torch.tensor(atom_mask, dtype=torch.bool).unsqueeze(0)  # [1, N]

        # Compute LDDT using existing infrastructure
        lddt_value = compute_lddt(
            all_atom_pred_pos=mobile_tensor,
            all_atom_positions=target_tensor,
            all_atom_mask=mask_tensor,
            per_residue=False,  # Return overall LDDT score
            **kwargs,
        )

        return float(lddt_value)

    def state_dict(self):
        """This state dict is optimized for storage, so it turns things to fp16 whenever
        possible and converts numpy arrays to lists for JSON serialization.
        """
        dct = {k: v for k, v in vars(self).items()}
        for k, v in dct.items():
            if isinstance(v, np.ndarray):
                match v.dtype:
                    case np.int64:
                        dct[k] = v.astype(np.int32).tolist()
                    case np.float64 | np.float32:
                        dct[k] = v.astype(np.float16).tolist()
                    case _:
                        dct[k] = v.tolist()
            elif isinstance(v, MolecularComplexMetadata):
                dct[k] = asdict(v)

        return dct

    def to_blob(self) -> bytes:
        return brotli.compress(msgpack.dumps(self.state_dict()), quality=5)

    @classmethod
    def from_state_dict(cls, dct):
        for k, v in dct.items():
            if isinstance(v, list) and k in [
                "atom_positions",
                "atom_elements",
                "token_to_atoms",
                "plddt",
            ]:
                dct[k] = np.array(v)

        for k, v in dct.items():
            if isinstance(v, np.ndarray):
                if k in ["atom_positions", "plddt"]:
                    dct[k] = v.astype(np.float32)
                elif k in ["token_to_atoms"]:
                    dct[k] = v.astype(np.int32)

        dct["metadata"] = MolecularComplexMetadata(**dct["metadata"])
        return cls(**dct)

    @classmethod
    def from_blob(cls, input: Path | str | io.BytesIO | bytes):
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
