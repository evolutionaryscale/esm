from typing import cast

import numpy as np

from examples.local_generate import get_sample_protein
from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    GenerationConfig,
)
from esm.sdk.forge import FoldForgeInferenceClient


def convert_none_to_nan(data):
    """Recursively convert None values in any deeply nested structure (e.g., list of lists of lists) to np.nan."""
    if isinstance(data, list):
        return [convert_none_to_nan(x) for x in data]
    elif data is None:
        return np.nan
    else:
        return data


def are_allclose_with_nan(A, B, rtol=1e-5, atol=1e-2):
    B = convert_none_to_nan(B)

    A = np.array(A)
    B = np.array(B)

    if A.shape != B.shape:
        raise ValueError("A and B must have the same shape")

    nan_mask_A = np.isnan(A)
    nan_mask_B = np.isnan(B)

    if not np.array_equal(nan_mask_A, nan_mask_B):
        return False

    return np.allclose(A[~nan_mask_A], B[~nan_mask_B], rtol=rtol, atol=atol)


def main(fold_client: FoldForgeInferenceClient, esm3_client: ESM3InferenceClient):
    # Folding
    protein = get_sample_protein()
    sequence_length = len(protein.sequence)  # type: ignore
    num_steps = int(sequence_length / 16)
    protein.coordinates = None
    protein.function_annotations = None
    protein.sasa = None
    # Folding with esm3 client
    folded_protein = cast(
        ESMProtein,
        esm3_client.generate(
            protein,
            GenerationConfig(
                track="structure", schedule="cosine", num_steps=num_steps, temperature=0
            ),
        ),
    )
    # Folding with folding client
    coordinates = fold_client.fold(
        "esm3",
        protein.sequence,  # type:ignore
        potential_sequence_of_concern=False,
    )
    assert are_allclose_with_nan(folded_protein.coordinates, coordinates)
