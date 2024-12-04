from esm.models.esmc import ESMC
from examples.local_generate import get_sample_protein
from esm.sdk.api import (
    ESMCInferenceClient,
    LogitsConfig,
    LogitsOutput,
)


def main(client: ESMCInferenceClient):
    # ================================================================
    # Example usage: one single protein
    # ================================================================
    protein = get_sample_protein()
    protein.coordinates = None
    protein.function_annotations = None
    protein.sasa = None

    # Use logits endpoint. Using bf16 for inference optimization
    protein_tensor = client.encode(protein)
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )
    assert isinstance(
        logits_output, LogitsOutput
    ), f"LogitsOutput was expected but got {logits_output}"
    assert (
        logits_output.logits is not None and logits_output.logits.sequence is not None
    )
    assert logits_output.embeddings is not None and logits_output.embeddings is not None


if __name__ == "__main__":
    main(ESMC.from_pretrained("esmc_300m"))
