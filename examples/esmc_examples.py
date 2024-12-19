from local_generate import get_sample_protein

from esm.models.esmc import ESMC
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
    output = client.logits(
        protein_tensor,
        LogitsConfig(sequence=True, return_embeddings=True, return_hidden_states=True),
    )
    assert isinstance(
        output, LogitsOutput
    ), f"LogitsOutput was expected but got {output}"
    assert output.logits is not None and output.logits.sequence is not None
    assert output.embeddings is not None
    assert output.hidden_states is not None
    print(
        f"Client returned logits with shape: {output.logits.sequence.shape}, embeddings with shape: {output.embeddings.shape}, and hidden states with shape {output.hidden_states.shape}"
    )


def raw_forward(model: ESMC):
    protein = get_sample_protein()
    assert protein.sequence is not None
    sequences = [protein.sequence, protein.sequence]

    # ================================================================
    # Example usage: directly use the model
    # ================================================================
    input_ids = model._tokenize(sequences)
    output = model(input_ids)
    logits, embeddings, hiddens = (
        output.sequence_logits,
        output.embeddings,
        output.hidden_states,
    )
    print(
        f"Raw model returned logits with shape: {logits.shape}, embeddings with shape: {embeddings.shape} and hidden states with shape {hiddens.shape}"
    )


if __name__ == "__main__":
    model = ESMC.from_pretrained("esmc_300m")
    main(model)
    raw_forward(model)
