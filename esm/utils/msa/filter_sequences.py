import os
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from esm.utils.system import run_subprocess_with_errorcheck


def greedy_select_indices(array, num_seqs: int, mode: str = "max") -> list[int]:
    """Greedily select sequences that either maximize or minimize hamming distance.

    Algorithm proposed in the MSA Transformer paper. Starting from the query sequence,
    iteratively add sequences to the list with the maximum (minimum) average Hamming
    distance to the existing set of sequences.

    Args:
        array (np.ndarray): Character array representing the sequences in the MSA
        num_seqs (int): Number of sequences to select.
        mode (str): Whether to maximize or minimize diversity. DO NOT pick 'min' unless
            you're doing it to prove a point for a paper.

    Returns:
        list[int]: List of indices to select from the array
    """
    assert mode in ("max", "min")
    depth = array.shape[0]
    if depth <= num_seqs:
        return list(range(depth))
    array = array.view(np.uint8)

    optfunc = np.argmax if mode == "max" else np.argmin
    all_indices = np.arange(depth)
    indices = [0]
    pairwise_distances = np.zeros((0, depth))
    for _ in range(num_seqs - 1):
        dist = cdist(array[indices[-1:]], array, "hamming")
        pairwise_distances = np.concatenate([pairwise_distances, dist])
        shifted_distance = np.delete(pairwise_distances, indices, axis=1).mean(0)
        shifted_index = optfunc(shifted_distance)
        index = np.delete(all_indices, indices)[shifted_index]
        indices.append(index)
    indices = sorted(indices)
    return indices


def hhfilter(
    sequences: list[str],
    seqid: int = 90,
    diff: int = 0,
    cov: int = 0,
    qid: int = 0,
    qsc: float = -20.0,
    binary: str = "hhfilter",
) -> list[int]:
    with tempfile.TemporaryDirectory(
        dir="/dev/shm" if os.path.exists("/dev/shm") else None
    ) as tempdirname:
        tempdir = Path(tempdirname)
        fasta_file = tempdir / "input.fasta"
        fasta_file.write_text(
            "\n".join(f">{i}\n{seq}" for i, seq in enumerate(sequences))
        )
        output_file = tempdir / "output.fasta"
        command = " ".join(
            [
                f"{binary}",
                f"-i {fasta_file}",
                "-M a3m",
                f"-o {output_file}",
                f"-id {seqid}",
                f"-diff {diff}",
                f"-cov {cov}",
                f"-qid {qid}",
                f"-qsc {qsc}",
            ]
        ).split(" ")
        run_subprocess_with_errorcheck(command, capture_output=True)
        with output_file.open() as f:
            indices = [int(line[1:].strip()) for line in f if line.startswith(">")]
        return indices
