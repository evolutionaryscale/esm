from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np


@dataclass
class Modification:
    position: int  # zero-indexed
    ccd: str


@dataclass
class ProteinInput:
    id: str | list[str]
    sequence: str
    modifications: list[Modification] | None = None


@dataclass
class RNAInput:
    id: str | list[str]
    sequence: str
    modifications: list[Modification] | None = None


@dataclass
class DNAInput:
    id: str | list[str]
    sequence: str
    modifications: list[Modification] | None = None


@dataclass
class LigandInput:
    id: str | list[str]
    smiles: str
    ccd: list[str] | None = None


@dataclass
class DistogramConditioning:
    chain_id: str
    distogram: np.ndarray


@dataclass
class PocketConditioning:
    binder_chain_id: str
    contacts: list[tuple[str, int]]


@dataclass
class CovalentBond:
    chain_id1: str
    res_idx1: int
    atom_idx1: int
    chain_id2: str
    res_idx2: int
    atom_idx2: int


@dataclass
class StructurePredictionInput:
    sequences: Sequence[ProteinInput | RNAInput | DNAInput | LigandInput]
    pocket: PocketConditioning | None = None
    distogram_conditioning: list[DistogramConditioning] | None = None
    covalent_bonds: list[CovalentBond] | None = None


def serialize_structure_prediction_input(all_atom_input: StructurePredictionInput):
    def create_chain_data(seq_input, chain_type: str) -> dict[str, Any]:
        chain_data: dict[str, Any] = {
            "sequence": seq_input.sequence,
            "id": seq_input.id,
            "type": chain_type,
        }
        if hasattr(seq_input, "modifications") and seq_input.modifications:
            mods = [
                {"position": mod.position, "ccd": mod.ccd}
                for mod in seq_input.modifications
            ]
            chain_data["modifications"] = mods
        return chain_data

    sequences = []
    for seq_input in all_atom_input.sequences:
        if isinstance(seq_input, ProteinInput):
            sequences.append(create_chain_data(seq_input, "protein"))
        elif isinstance(seq_input, RNAInput):
            sequences.append(create_chain_data(seq_input, "rna"))
        elif isinstance(seq_input, DNAInput):
            sequences.append(create_chain_data(seq_input, "dna"))
        elif isinstance(seq_input, LigandInput):
            sequences.append(
                {
                    "smiles": seq_input.smiles,
                    "id": seq_input.id,
                    "ccd": seq_input.ccd,
                    "type": "ligand",
                }
            )
        else:
            raise ValueError(f"Unsupported sequence input type: {type(seq_input)}")

    result: dict[str, Any] = {"sequences": sequences}

    # Add covalent bonds if present
    if all_atom_input.covalent_bonds is not None:
        result["covalent_bonds"] = [
            {
                "chain_id1": bond.chain_id1,
                "res_idx1": bond.res_idx1,
                "atom_idx1": bond.atom_idx1,
                "chain_id2": bond.chain_id2,
                "res_idx2": bond.res_idx2,
                "atom_idx2": bond.atom_idx2,
            }
            for bond in all_atom_input.covalent_bonds
        ]

    return result
