"""Device dispatch helpers for CUDA, XPU, and CPU."""

from __future__ import annotations

import torch

_ACCELERATED_AUTOCAST_TYPES = ("cuda", "xpu")


def resolve_default_device() -> torch.device:
    """Best available accelerator, else CPU."""
    if torch.accelerator.is_available():
        device = torch.accelerator.current_accelerator()
        if device is None:
            # is_available() said yes, current_accelerator() said no -- that's
            # a torch bug or a broken install, not a case to silently fall
            # back to CPU on.
            raise RuntimeError(
                "torch.accelerator.is_available() is True but "
                "current_accelerator() returned None."
            )
        return device
    return torch.device("cpu")


def supports_amp_autocast(device_type: str) -> bool:
    """Whether device_type supports bf16/fp16 autocast (cuda, xpu do; mps/cpu don't the same way)."""
    return device_type in _ACCELERATED_AUTOCAST_TYPES


def empty_cache(device_type: str) -> None:
    if device_type == "cuda":
        torch.cuda.empty_cache()
    elif device_type == "xpu":
        torch.xpu.empty_cache()


def synchronize(device_type: str) -> None:
    if device_type == "cuda":
        torch.cuda.synchronize()
    elif device_type == "xpu":
        torch.xpu.synchronize()


def manual_seed_all(device_type: str, seed: int) -> None:
    if device_type == "cuda":
        torch.cuda.manual_seed_all(seed)
    elif device_type == "xpu":
        torch.xpu.manual_seed_all(seed)


def get_rng_state_all(device_type: str) -> list[torch.Tensor] | None:
    if device_type == "cuda":
        return torch.cuda.get_rng_state_all()
    elif device_type == "xpu":
        return torch.xpu.get_rng_state_all()
    return None


def set_rng_state_all(device_type: str, states: list[torch.Tensor]) -> None:
    if device_type == "cuda":
        torch.cuda.set_rng_state_all(states)
    elif device_type == "xpu":
        torch.xpu.set_rng_state_all(states)
