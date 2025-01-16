import os
from collections import defaultdict
from io import BytesIO
from typing import Any, ContextManager, Sequence, TypeVar
from warnings import warn

import huggingface_hub
import numpy as np
import torch
import zstd

from esm.utils.constants.esm3 import CHAIN_BREAK_STR
from esm.utils.types import FunctionAnnotation

MAX_SUPPORTED_DISTANCE = 1e6


TSequence = TypeVar("TSequence", bound=Sequence)


def slice_python_object_as_numpy(
    obj: TSequence, idx: int | list[int] | slice | np.ndarray
) -> TSequence:
    """
    Slice a python object (like a list, string, or tuple) as if it was a numpy object.

    Example:
        >>> obj = "ABCDE"
        >>> slice_python_object_as_numpy(obj, [1, 3, 4])
        "BDE"

        >>> obj = [1, 2, 3, 4, 5]
        >>> slice_python_object_as_numpy(obj, np.arange(5) < 3)
        [1, 2, 3]
    """
    if isinstance(idx, int):
        idx = [idx]

    if isinstance(idx, np.ndarray) and idx.dtype == bool:
        sliced_obj = [obj[i] for i in np.where(idx)[0]]
    elif isinstance(idx, slice):
        sliced_obj = obj[idx]
    else:
        sliced_obj = [obj[i] for i in idx]

    match obj, sliced_obj:
        case str(), list():
            sliced_obj = "".join(sliced_obj)
        case _:
            sliced_obj = obj.__class__(sliced_obj)  # type: ignore

    return sliced_obj  # type: ignore


def rbf(values, v_min, v_max, n_bins=16):
    """
    Returns RBF encodings in a new dimension at the end.
    """
    rbf_centers = torch.linspace(
        v_min, v_max, n_bins, device=values.device, dtype=values.dtype
    )
    rbf_centers = rbf_centers.view([1] * len(values.shape) + [-1])
    rbf_std = (v_max - v_min) / n_bins
    z = (values.unsqueeze(-1) - rbf_centers) / rbf_std
    return torch.exp(-(z**2))


def batched_gather(data, inds, dim=0, no_batch_dims=0):
    ranges = []
    for i, s in enumerate(data.shape[:no_batch_dims]):
        r = torch.arange(s)
        r = r.view(*(*((1,) * i), -1, *((1,) * (len(inds.shape) - i - 1))))
        ranges.append(r)

    remaining_dims = [slice(None) for _ in range(len(data.shape) - no_batch_dims)]
    remaining_dims[dim - no_batch_dims if dim >= 0 else dim] = inds
    ranges.extend(remaining_dims)
    return data[ranges]


def node_gather(s: torch.Tensor, edges: torch.Tensor) -> torch.Tensor:
    return batched_gather(s.unsqueeze(-3), edges, -2, no_batch_dims=len(s.shape) - 1)


def knn_graph(
    coords: torch.Tensor,
    coord_mask: torch.Tensor,
    padding_mask: torch.Tensor,
    sequence_id: torch.Tensor,
    *,
    no_knn: int,
):
    L = coords.shape[-2]
    num_by_dist = min(no_knn, L)
    device = coords.device

    coords = coords.nan_to_num()
    coord_mask = ~(coord_mask[..., None, :] & coord_mask[..., :, None])
    padding_pairwise_mask = padding_mask[..., None, :] | padding_mask[..., :, None]
    if sequence_id is not None:
        padding_pairwise_mask |= torch.unsqueeze(sequence_id, 1) != torch.unsqueeze(
            sequence_id, 2
        )
    dists = (coords.unsqueeze(-2) - coords.unsqueeze(-3)).norm(dim=-1)
    arange = torch.arange(L, device=device)
    seq_dists = (arange.unsqueeze(-1) - arange.unsqueeze(-2)).abs()
    # We only support up to a certain distance, above that, we use sequence distance
    # instead. This is so that when a large portion of the structure is masked out,
    # the edges are built according to sequence distance.
    max_dist = MAX_SUPPORTED_DISTANCE
    torch._assert_async((dists[~coord_mask] < max_dist).all())
    struct_then_seq_dist = (
        seq_dists.to(dists.dtype)
        .mul(1e2)
        .add(max_dist)
        .where(coord_mask, dists)
        .masked_fill(padding_pairwise_mask, torch.inf)
    )
    dists, edges = struct_then_seq_dist.sort(dim=-1, descending=False)
    # This is a L x L tensor, where we index by rows first,
    # and columns are the edges we should pick.
    chosen_edges = edges[..., :num_by_dist]
    chosen_mask = dists[..., :num_by_dist].isfinite()
    return chosen_edges, chosen_mask


def stack_variable_length_tensors(
    sequences: Sequence[torch.Tensor],
    constant_value: int | float = 0,
    dtype: torch.dtype | None = None,
) -> torch.Tensor:
    """Automatically stack tensors together, padding variable lengths with the
    value in constant_value. Handles an arbitrary number of dimensions.

    Examples:
        >>> tensor1, tensor2 = torch.ones([2]), torch.ones([5])
        >>> stack_variable_length_tensors(tensor1, tensor2)
        tensor of shape [2, 5]. First row is [1, 1, 0, 0, 0]. Second row is all ones.

        >>> tensor1, tensor2 = torch.ones([2, 4]), torch.ones([5, 3])
        >>> stack_variable_length_tensors(tensor1, tensor2)
        tensor of shape [2, 5, 4]
    """
    batch_size = len(sequences)
    shape = [batch_size] + np.max([seq.shape for seq in sequences], 0).tolist()

    if dtype is None:
        dtype = sequences[0].dtype
    device = sequences[0].device

    array = torch.full(shape, constant_value, dtype=dtype, device=device)
    for arr, seq in zip(array, sequences):
        arrslice = tuple(slice(dim) for dim in seq.shape)
        arr[arrslice] = seq

    return array


def unbinpack(
    tensor: torch.Tensor, sequence_id: torch.Tensor | None, pad_value: int | float
):
    """
    Args:
        tensor (Tensor): [B, L, ...]

    Returns:
        Tensor: [B_unbinpacked, L_unbinpack, ...]
    """
    if sequence_id is None:
        return tensor

    unpacked_tensors = []
    num_sequences = sequence_id.max(dim=-1).values + 1
    for batch_idx, (batch_seqid, batch_num_sequences) in enumerate(
        zip(sequence_id, num_sequences)
    ):
        for seqid in range(batch_num_sequences):
            mask = batch_seqid == seqid
            unpacked = tensor[batch_idx, mask]
            unpacked_tensors.append(unpacked)
    return stack_variable_length_tensors(unpacked_tensors, pad_value)


def fp32_autocast_context(device_type: str) -> ContextManager[torch.amp.autocast]:  # type: ignore
    """
    Returns an autocast context manager that disables downcasting by AMP.

    Args:
        device_type: The device type ('cpu' or 'cuda')

    Returns:
        An autocast context manager with the specified behavior.
    """
    if device_type == "cpu":
        return torch.amp.autocast(device_type, enabled=False)  # type: ignore
    elif device_type == "cuda":
        return torch.amp.autocast(device_type, dtype=torch.float32)  # type: ignore
    else:
        raise ValueError(f"Unsupported device type: {device_type}")


def merge_ranges(ranges: list[range], merge_gap_max: int | None = None) -> list[range]:
    """Merge overlapping ranges into sorted, non-overlapping segments.

    Args:
        ranges: collection of ranges to merge.
        merge_gap_max: optionally merge neighboring ranges that are separated by a gap
          no larger than this size.
    Returns:
        non-overlapping ranges merged from the inputs, sorted by position.
    """
    ranges = sorted(ranges, key=lambda r: r.start)
    merge_gap_max = merge_gap_max if merge_gap_max is not None else 0
    assert merge_gap_max >= 0, f"Invalid merge_gap_max: {merge_gap_max}"

    merged = []
    for r in ranges:
        if not merged:
            merged.append(r)
        else:
            last = merged[-1]
            if last.stop + merge_gap_max >= r.start:
                merged[-1] = range(last.start, max(last.stop, r.stop))
            else:
                merged.append(r)
    return merged


def merge_annotations(
    annotations: list[FunctionAnnotation], merge_gap_max: int | None = None
) -> list[FunctionAnnotation]:
    """Merges annotations into non-overlapping segments.

    Args:
        annotations: annotations to merge.
        merge_gap_max: optionally merge neighboring ranges that are separated by a gap
          no larger than this size.
    Returns:
        non-overlapping annotations with gaps merged.
    """
    grouped: dict[str, list[range]] = defaultdict(list)
    for a in annotations:
        # +1 since FunctionAnnotation.end is inlcusive.
        grouped[a.label].append(range(a.start, a.end + 1))

    merged = []
    for label, ranges in grouped.items():
        merged_ranges = merge_ranges(ranges, merge_gap_max=merge_gap_max)
        for range_ in merged_ranges:
            annotation = FunctionAnnotation(
                label=label,
                start=range_.start,
                end=range_.stop - 1,  # convert range.stop exclusive -> inclusive.
            )
            merged.append(annotation)
    return merged


def replace_inf(data):
    if data is None:
        return None
    array = np.array(data, dtype=np.float32, copy=False)
    array = np.where(np.isinf(array), -1, array)
    return array.tolist()


def maybe_tensor(x, convert_none_to_nan: bool = False) -> torch.Tensor | None:
    if x is None:
        return None
    if convert_none_to_nan:
        x = np.array(x, copy=False, dtype=np.float32)
        x = np.where(x is None, np.nan, x)
    return torch.tensor(x)


def maybe_list(x, convert_nan_to_none: bool = False) -> list | None:
    if x is None:
        return None
    if not convert_nan_to_none:
        return x.tolist()
    nan_mask = torch.isnan(x)
    np_arr = x.cpu().numpy().astype(object)
    np_arr[nan_mask.cpu().numpy()] = None
    return np_arr.tolist()


def huggingfacehub_login():
    """Authenticates with the Hugging Face Hub using the HF_TOKEN environment
    variable, else by prompting the user"""
    token = os.environ.get("HF_TOKEN")
    huggingface_hub.login(token=token)


def get_chainbreak_boundaries_from_sequence(sequence: Sequence[str]) -> np.ndarray:
    chain_boundaries = [0]
    for i, aa in enumerate(sequence):
        if aa == CHAIN_BREAK_STR:
            if i == (len(sequence) - 1):
                raise ValueError(
                    "Encountered chain break token at end of sequence, this is unexpected."
                )
            if i == (len(sequence) - 2):
                warn(
                    "Encountered chain break token at penultimate position, this is unexpected."
                )
            chain_boundaries.append(i)
            chain_boundaries.append(i + 1)
    chain_boundaries.append(len(sequence))
    assert len(chain_boundaries) % 2 == 0
    chain_boundaries = np.array(chain_boundaries).reshape(-1, 2)
    return chain_boundaries


def deserialize_tensors(b: bytes) -> Any:
    buf = BytesIO(zstd.ZSTD_uncompress(b))
    d = torch.load(buf, map_location="cpu", weights_only=False)
    return d
