"""Inference-only fp32 RoPE kernel for ESMC on Intel XPU."""

import torch
import triton
import triton.language as tl


@triton.jit
def _rotary_kernel(
    x_ptr,
    cos_ptr,
    sin_ptr,
    out_ptr,
    B: tl.constexpr,
    S: tl.constexpr,
    H: tl.constexpr,
    D: tl.constexpr,
    RO: tl.constexpr,
    X_S0: tl.constexpr,
    X_S1: tl.constexpr,
    X_S2: tl.constexpr,
    X_S3: tl.constexpr,
    C_S0: tl.constexpr,
    C_S1: tl.constexpr,
    S_S0: tl.constexpr,
    S_S1: tl.constexpr,
    BLOCK: tl.constexpr,
):
    index = tl.program_id(0) * BLOCK + tl.arange(0, BLOCK)
    total = B * S * H * D
    mask = index < total

    d = index % D
    h = (index // D) % H
    s = (index // (D * H)) % S
    b = index // (D * H * S)
    x_offset = b * X_S0 + s * X_S1 + h * X_S2 + d * X_S3
    value = tl.load(x_ptr + x_offset, mask=mask, other=0.0).to(tl.float32)

    rotate_mask = mask & (d < RO)
    frequency = tl.where(d < RO // 2, d, d - RO // 2)
    cos_value = tl.load(
        cos_ptr + s * C_S0 + frequency * C_S1, mask=rotate_mask, other=0.0
    ).to(tl.float32)
    sin_value = tl.load(
        sin_ptr + s * S_S0 + frequency * S_S1, mask=rotate_mask, other=0.0
    ).to(tl.float32)
    pair_d = tl.where(d < RO // 2, d + RO // 2, d - RO // 2)
    pair_offset = b * X_S0 + s * X_S1 + h * X_S2 + pair_d * X_S3
    pair_value = tl.load(x_ptr + pair_offset, mask=rotate_mask, other=0.0).to(
        tl.float32
    )
    pair_sign = tl.where(d < RO // 2, -1.0, 1.0)
    rotated = value * cos_value + (pair_value * pair_sign) * sin_value
    result = tl.where(d < RO, rotated, value)
    tl.store(out_ptr + index, result, mask=mask)


def apply_rotary_emb_xpu(
    x: torch.Tensor, cos: torch.Tensor, sin: torch.Tensor
) -> torch.Tensor:
    """Apply ESMC's non-interleaved RoPE in fp32 inference mode."""
    if torch.is_grad_enabled():
        raise RuntimeError("XPU Triton RoPE is only supported with gradients disabled")
    if x.device.type != "xpu" or cos.device != x.device or sin.device != x.device:
        raise ValueError("XPU Triton RoPE requires all inputs on the same XPU")
    if x.dtype != torch.float32 or cos.dtype != x.dtype or sin.dtype != x.dtype:
        raise ValueError("XPU Triton RoPE requires fp32 inputs")
    if x.ndim != 4 or cos.ndim != 2 or sin.shape != cos.shape:
        raise ValueError("unexpected RoPE input shapes")

    batch, sequence, heads, dim = x.shape
    rotary_dim = cos.shape[-1] * 2
    if rotary_dim > dim or cos.shape[0] < sequence:
        raise ValueError("RoPE dimensions exceed the input shape")

    output = torch.empty(
        x.shape, device=x.device, dtype=x.dtype, memory_format=torch.contiguous_format
    )
    block = 128
    _rotary_kernel[(triton.cdiv(x.numel(), block),)](
        x,
        cos,
        sin,
        output,
        batch,  # ty:ignore[invalid-argument-type]
        sequence,  # ty:ignore[invalid-argument-type]
        heads,  # ty:ignore[invalid-argument-type]
        dim,  # ty:ignore[invalid-argument-type]
        rotary_dim,  # ty:ignore[invalid-argument-type]
        x.stride(0),  # ty:ignore[invalid-argument-type]
        x.stride(1),  # ty:ignore[invalid-argument-type]
        x.stride(2),  # ty:ignore[invalid-argument-type]
        x.stride(3),  # ty:ignore[invalid-argument-type]
        cos.stride(0),  # ty:ignore[invalid-argument-type]
        cos.stride(1),  # ty:ignore[invalid-argument-type]
        sin.stride(0),  # ty:ignore[invalid-argument-type]
        sin.stride(1),  # ty:ignore[invalid-argument-type]
        BLOCK=block,  # ty:ignore[invalid-argument-type]
        num_warps=4,  # ty:ignore[unknown-argument]
    )
    return output
