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
class StructurePredictionInput:
    sequences: Sequence[ProteinInput | RNAInput | DNAInput | LigandInput]
    pocket: PocketConditioning | None = None
    distogram_conditioning: list[DistogramConditioning] | None = None


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

    return {"sequences": sequences}
