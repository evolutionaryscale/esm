from typing import Callable

import torch
import torch.nn as nn

from esm.models.esm3 import ESM3
from esm.models.function_decoder import FunctionTokenDecoder
from esm.models.vqvae import (
    StructureTokenDecoder,
    StructureTokenEncoder,
)
from esm.utils.constants.esm3 import data_root
from esm.utils.constants.models import (
    ESM3_FUNCTION_DECODER_V0,
    ESM3_OPEN_SMALL,
    ESM3_STRUCTURE_DECODER_V0,
    ESM3_STRUCTURE_ENCODER_V0,
)

ModelBuilder = Callable[[torch.device | str], nn.Module]


def ESM3_sm_open_v0(device: torch.device | str = "cpu"):
    model = (
        ESM3(
            d_model=1536,
            n_heads=24,
            v_heads=256,
            n_layers=48,
            structure_encoder_name=ESM3_STRUCTURE_ENCODER_V0,
            structure_decoder_name=ESM3_STRUCTURE_DECODER_V0,
            function_decoder_name=ESM3_FUNCTION_DECODER_V0,
        )
        .to(device)
        .eval()
    )
    state_dict = torch.load(
        data_root() / "data/weights/esm3_sm_open_v1.pth", map_location=device
    )
    model.load_state_dict(state_dict)
    return model


def ESM3_structure_encoder_v0(device: torch.device | str = "cpu"):
    model = (
        StructureTokenEncoder(
            d_model=1024, n_heads=1, v_heads=128, n_layers=2, d_out=128, n_codes=4096
        )
        .to(device)
        .eval()
    )
    state_dict = torch.load(
        data_root() / "data/weights/esm3_structure_encoder_v0.pth", map_location=device
    )
    model.load_state_dict(state_dict)
    return model


def ESM3_structure_decoder_v0(device: torch.device | str = "cpu"):
    model = (
        StructureTokenDecoder(d_model=1280, n_heads=20, n_layers=30).to(device).eval()
    )
    state_dict = torch.load(
        data_root() / "data/weights/esm3_structure_decoder_v0.pth", map_location=device
    )
    model.load_state_dict(state_dict)
    return model


def ESM3_function_decoder_v0(device: torch.device | str = "cpu"):
    model = FunctionTokenDecoder().to(device).eval()
    state_dict = torch.load(
        data_root() / "data/weights/esm3_function_decoder_v0.pth", map_location=device
    )
    model.load_state_dict(state_dict)
    return model


LOCAL_MODEL_REGISTRY: dict[str, ModelBuilder] = {
    ESM3_OPEN_SMALL: ESM3_sm_open_v0,
    ESM3_STRUCTURE_ENCODER_V0: ESM3_structure_encoder_v0,
    ESM3_STRUCTURE_DECODER_V0: ESM3_structure_decoder_v0,
    ESM3_FUNCTION_DECODER_V0: ESM3_function_decoder_v0,
}


def load_local_model(model_name: str, device: torch.device | str = "cpu") -> nn.Module:
    if model_name not in LOCAL_MODEL_REGISTRY:
        raise ValueError(f"Model {model_name} not found in local model registry.")
    return LOCAL_MODEL_REGISTRY[model_name](device)


# Register custom versions of ESM3 for use with the local inference API
def register_local_model(model_name: str, model_builder: ModelBuilder) -> None:
    LOCAL_MODEL_REGISTRY[model_name] = model_builder
