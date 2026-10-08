"""Synthetic ESM3 XPU tests: small random-init model, no real checkpoint.

Mirrors tests/oss_pytests/test_output_attentions.py's make_esm3() pattern,
checking that ESM3's forward pass dispatches correctly on XPU and agrees
with the CPU reference.
"""

from unittest.mock import MagicMock

import pytest
import torch

from esm import pretrained
from esm.models.esm3 import ESM3, ESMOutput
from esm.utils.constants import esm3 as C

B, L = 2, 10
D_MODEL = 32
N_HEADS = 4
N_LAYERS = 2


def make_esm3() -> ESM3:
    tokenizers = MagicMock()
    tokenizers.sequence.mask_token_id = 32

    return ESM3(
        d_model=D_MODEL,
        n_heads=N_HEADS,
        v_heads=4,
        n_layers=N_LAYERS,
        structure_encoder_fn=MagicMock(),
        structure_decoder_fn=MagicMock(),
        function_decoder_fn=MagicMock(),
        tokenizers=tokenizers,
    )


def test_from_pretrained_accepts_xpu_device_string(monkeypatch):
    model = make_esm3()
    loaded_devices: list[torch.device] = []

    def load_local_model(model_name: str, device: torch.device) -> ESM3:
        loaded_devices.append(device)
        return model

    monkeypatch.setattr(pretrained, "load_local_model", load_local_model)
    assert ESM3.from_pretrained(device="xpu") is model
    assert loaded_devices == [torch.device("xpu")]


@pytest.mark.xpu
def test_cpu_and_xpu_agree():
    """Same weights, same input, both devices, fp32 on each side."""
    torch.manual_seed(0)
    model = make_esm3().eval()
    sequence_tokens = torch.full((B, L), C.SEQUENCE_MASK_TOKEN, dtype=torch.long)
    sequence_tokens[:, 0] = C.SEQUENCE_BOS_TOKEN
    sequence_tokens[:, -1] = C.SEQUENCE_EOS_TOKEN

    with torch.no_grad():
        expected: ESMOutput = model(sequence_tokens=sequence_tokens)

        model_xpu = model.to("xpu")
        actual: ESMOutput = model_xpu(sequence_tokens=sequence_tokens.to("xpu"))

    torch.testing.assert_close(
        actual.sequence_logits.cpu(), expected.sequence_logits, atol=1e-5, rtol=1e-5
    )
