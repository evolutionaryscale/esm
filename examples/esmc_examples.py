from esm.models.esmc import ESMC
from esm.sdk.api import ESMCInferenceClient, ESMProtein, LogitsConfig, LogitsOutput


def main(client: ESMCInferenceClient):
    # ================================================================
    # Example usage: one single protein
    # ================================================================
    protein = ESMProtein(sequence="AAAAA")

    # Use logits endpoint. Using bf16 for inference optimization
    protein_tensor = client.encode(protein)
    output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )
    assert isinstance(
        output, LogitsOutput
    ), f"LogitsOutput was expected but got {output}"
    assert output.logits is not None and output.logits.sequence is not None
    assert output.embeddings is not None and output.embeddings is not None
    print(
        f"Client returned logits with shape: {output.logits.sequence.shape} and embeddings with shape: {output.embeddings.shape}"
    )


def raw_forward(model: ESMC):
    protein = ESMProtein(sequence="AAAAA")
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
