import pytest
import torch

from esm.sdk.api import (
    ESMProtein,
    ESMProteinTensor,
    GenerationConfig,
)
from evolutionaryscale.utils.env import ModelName
from evolutionaryscale.utils.remote_inference.api_v1 import (
    ESM3RemoteModelInferenceClient,
)
from projects.forge.fastapi.utils.model import _load_esm_model


@pytest.fixture()
def esm3_remote_inference_client():
    model = _load_esm_model(
        ModelName.ESM3_TINY_DEV, distributed_model=False, load_function_decoder=False
    )
    client = ESM3RemoteModelInferenceClient(
        model,
        tokenizers=model.tokenizers,
        device=torch.device("cuda"),
        enable_batched_runner=False,
    )
    return client


@pytest.mark.gpu
def test_chain_break_tokens(esm3_remote_inference_client):
    tokenizer = esm3_remote_inference_client.tokenizers.sequence
    # 3 separate chains with 2 chainbreak tokens.
    sequence_with_chain_breaks = torch.tensor(
        [
            tokenizer.bos_token_id,
            20,
            20,
            20,
            20,
            tokenizer.chain_break_token_id,
            21,
            21,
            21,
            tokenizer.chain_break_token_id,
            22,
            22,
            22,
            tokenizer.eos_token_id,
        ]
    )
    protein = esm3_remote_inference_client.generate(
        ESMProteinTensor(sequence=sequence_with_chain_breaks),
        # There are 10 tokens that actually need to be sampled.
        GenerationConfig(track="structure", num_steps=10),
    )

    assert isinstance(protein, ESMProteinTensor)
    assert protein.structure is not None


@pytest.mark.gpu
def test_num_decoding_steps_more_than_mask_tokens(esm3_remote_inference_client):
    protein = esm3_remote_inference_client.generate(
        esm3_remote_inference_client.encode(
            ESMProtein(sequence="CDEFG")
        ),  # sequence of 5.
        GenerationConfig(track="structure", num_steps=10),  # use 10 decoding steps.
    )
    # Client should handle over-specification of decoding steps.
    # TODO: This should be a warning.
    assert isinstance(protein, ESMProteinTensor)
    assert protein.structure is not None


@pytest.mark.gpu
def test_num_decoding_steps_more_than_mask_tokens_batched(esm3_remote_inference_client):
    protein_list = esm3_remote_inference_client.batch_generate(
        inputs=[
            esm3_remote_inference_client.encode(ESMProtein(sequence="CDEFG")),
            esm3_remote_inference_client.encode(ESMProtein(sequence="ABCDEFG")),
            esm3_remote_inference_client.encode(ESMProtein(sequence="AB__EFG")),
        ],
        configs=[
            GenerationConfig(track="structure", num_steps=10),
            GenerationConfig(track="structure", num_steps=3),
            GenerationConfig(track="sequence", num_steps=20),
        ],
    )
    # Client should handle over-specification of decoding steps.
    # TODO: This should be a warning.
    assert isinstance(protein_list[0], ESMProteinTensor)
    assert protein_list[0].structure is not None
    assert isinstance(protein_list[1], ESMProteinTensor)
    assert protein_list[1].structure is not None
    assert isinstance(protein_list[2], ESMProteinTensor)
    assert protein_list[2].sequence is not None


@pytest.mark.gpu
def test_encode_chainbreak_token(esm3_remote_inference_client):
    protein = esm3_remote_inference_client.encode(ESMProtein(sequence="MSTNP|KPQKK"))
    assert isinstance(protein, ESMProteinTensor)
    assert protein.sequence is not None
    assert (
        protein.sequence[6]
        == esm3_remote_inference_client.tokenizers.sequence.chain_break_token_id
    )


@pytest.mark.gpu
def test_generation_with_chainbreak_token(esm3_remote_inference_client):
    chainbreak_sequence = torch.tensor(
        [
            esm3_remote_inference_client.tokenizers.sequence.bos_token_id,
            20,
            8,
            11,
            17,
            14,
            esm3_remote_inference_client.tokenizers.sequence.chain_break_token_id,
            15,
            14,
            16,
            15,
            15,
            esm3_remote_inference_client.tokenizers.sequence.eos_token_id,
        ]
    )

    protein = esm3_remote_inference_client.generate(
        ESMProteinTensor(sequence=chainbreak_sequence),
        GenerationConfig(track="structure", num_steps=1),
    )
    # Can't specify more decoding steps than masks available.
    assert isinstance(protein, ESMProteinTensor)
    assert protein.structure is not None
    assert (
        protein.structure[6]
        == esm3_remote_inference_client.tokenizers.structure.chain_break_token_id
    )
