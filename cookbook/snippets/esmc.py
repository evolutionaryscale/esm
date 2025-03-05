import os

from esm.models.esmc import ESMC
from esm.sdk import client
from esm.sdk.api import (
    ESMCInferenceClient,
    ESMProtein,
    ESMProteinTensor,
    LogitsConfig,
    LogitsOutput,
)
from esm.sdk.forge import ESM3ForgeInferenceClient


def main(client: ESMCInferenceClient | ESM3ForgeInferenceClient):
    # ================================================================
    # Example usage: one single protein
    # ================================================================
    protein = ESMProtein(sequence="AAAAA")

    # Use logits endpoint. Using bf16 for inference optimization
    protein_tensor = client.encode(protein)
    assert isinstance(
        protein_tensor, ESMProteinTensor
    ), f"Expected ESMProteinTensor but got error: {protein_tensor}"
    output = client.logits(
        protein_tensor,
        LogitsConfig(sequence=True, return_embeddings=True, return_hidden_states=True),
    )
    assert isinstance(
        output, LogitsOutput
    ), f"LogitsOutput was expected but got error: {output}"
    assert output.logits is not None and output.logits.sequence is not None
    assert output.embeddings is not None
    assert output.hidden_states is not None
    print(
        f"Client returned logits with shape: {output.logits.sequence.shape}, embeddings with shape: {output.embeddings.shape}, and hidden states with shape {output.hidden_states.shape}"
    )

    # request a specific hidden layer.
    assert isinstance(
        protein_tensor, ESMProteinTensor
    ), f"Expected ESMProteinTensor but got error: {protein_tensor}"
    output = client.logits(
        protein_tensor, LogitsConfig(return_hidden_states=True, ith_hidden_layer=1)
    )
    assert isinstance(
        output, LogitsOutput
    ), f"LogitsOutput was expected but got error: {output}"
    assert output.hidden_states is not None
    print(f"Client returned hidden states with shape {output.hidden_states.shape}")


def raw_forward(model: ESMC):
    protein = ESMProtein(sequence="AAAAA")
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
    if os.environ.get("ESM_API_KEY", ""):
        print("ESM_API_KEY found. Trying to use model from Forge...")
        main(client(model="esmc-300m-2024-12"))
    else:
        print("No ESM_API_KEY found. Trying to load model locally...")
        print(
            "TO try this script with a Forge API, please run ESM_API_KEY=your_api_key python esm3.py"
        )
        main(ESMC.from_pretrained("esm3_sm_open_v1"))
        model = ESMC.from_pretrained("esmc_300m")
        main(model)
        raw_forward(model)
