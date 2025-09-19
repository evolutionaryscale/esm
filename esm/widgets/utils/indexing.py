import numpy as np

from esm.utils.structure.protein_chain import ProteinChain

ZERO_INDEX = "Zero index"
PDB_INDEX = "PDB index"

PDB_INDEX_SUFFIX = "[PDB Index]"


def get_pdb_index_min_max(protein_chain: ProteinChain) -> tuple[int, int]:
    residue_index = protein_chain.residue_index
    valid_residue_index = residue_index[residue_index != -1]
    return min(valid_residue_index), max(valid_residue_index)


def pdb_index_to_zero_index(residue_index: int, protein_chain: ProteinChain) -> int:
    # Find the first position equal to residue_index
    pos = np.argwhere(residue_index == protein_chain.residue_index)
    if len(pos) == 0:
        raise ValueError(f"Residue index {residue_index} not found in protein chain")
    return pos[0][0]


def zero_index_to_pdb_index(zero_index: int, protein_chain: ProteinChain) -> int:
    return protein_chain.residue_index[zero_index]


def zero_range_to_pdb_range(
    zero_range: tuple[int, int], protein_chain: ProteinChain
) -> tuple[int, int]:
    return (
        zero_index_to_pdb_index(zero_range[0], protein_chain),
        zero_index_to_pdb_index(zero_range[1], protein_chain),
    )


def pdb_range_to_zero_range(
    pdb_range: tuple[int, int], protein_chain: ProteinChain
) -> tuple[int, int]:
    return (
        pdb_index_to_zero_index(pdb_range[0], protein_chain),
        pdb_index_to_zero_index(pdb_range[1], protein_chain),
    )
