import attr
import torch
import torch.nn.functional as F

from esm.sdk.api import (
    SamplingConfig,
    SamplingTrackConfig,
)
from esm.tokenization import (
    TokenizerCollection,
    get_invalid_tokenizer_ids,
)
from esm.tokenization.function_tokenizer import (
    InterProQuantizedTokenizer,
)
from esm.utils.constants.esm3 import MAX_RESIDUE_ANNOTATIONS


def get_default_sampling_config(tokenizers: TokenizerCollection) -> SamplingConfig:
    tracks = [f.name for f in attr.fields(SamplingConfig)]
    sampling_config = SamplingConfig()
    for current_track in tracks:
        setattr(
            sampling_config,
            current_track,
            SamplingTrackConfig(
                invalid_ids=get_invalid_tokenizer_ids(
                    getattr(tokenizers, current_track)
                ),
                temperature=1.0,
                top_p=1.0,
                # TODO: Add different mask and padding tokens for all tracks
                # Some tracks have the same pad and mask, which causes ambiguity when sampling
                only_sample_masked_tokens=current_track
                not in ["secondary_structure", "sasa", "function"],
            ),
        )
    return sampling_config


def sample_logits(
    logits: torch.Tensor,
    temperature: float | torch.Tensor,
    top_p: float | torch.Tensor = 1.0,
):
    """Default sampling from logits.

    Args:
        logits is shape (..., vocab_size)
        temperature is broadcastable to (...)
    """

    if top_p < 1.0:
        logits = top_p_logits(logits, top_p=top_p)

    temperature = _tensorize_like(temperature, logits)

    if torch.all(temperature == 0):
        ids = logits.argmax(-1)
        return ids

    assert not torch.any(temperature == 0), "Partial temperature 0 not supported."

    batch_dims = logits.size()[:-1]
    logits = logits.reshape(-1, logits.shape[-1])

    # Sample from all logits
    probs = F.softmax(logits / temperature[..., None], dim=-1)
    ids = torch.multinomial(probs, 1).squeeze(1)

    ids = ids.reshape(*batch_dims)
    return ids


def sample_function_logits(
    logits: torch.Tensor,
    tokenizer: InterProQuantizedTokenizer,
    top_p: float | torch.Tensor = 1.0,
    temperature: float | torch.Tensor = 1.0,
    p_none_threshold: float = 0.05,
) -> tuple[torch.Tensor, torch.Tensor]:
    [L, D, V] = logits.shape
    assert D == tokenizer.depth

    if top_p < 1.0:
        logits = top_p_logits(logits, top_p=top_p)

    temperature = torch.ones_like(logits[..., 0]) * temperature

    log_p = F.log_softmax(logits / temperature[..., None], dim=-1)  # (L, D, V)

    # Choose which positions have no predicted function.
    log_p_nones = log_p[..., tokenizer.vocab_to_index["<none>"]]  # (L, D)
    p_none = torch.exp(log_p_nones).mean(dim=-1)  # "Ensemble of <none> predictions"
    where_none = p_none > p_none_threshold  # (L, )

    # Set probability of <none> to 0 for all not-none positions
    none_index = tokenizer.vocab_to_index["<none>"]
    log_p[~where_none, :, none_index] = -torch.inf

    ids = torch.argmax(log_p, dim=-1)  # (L, D)
    ids[where_none, :] = tokenizer.vocab_to_index["<none>"]

    return ids, log_p


def sample_residue_annotation_logits(
    logits: torch.Tensor, annotation_threshold: float = 0.5
) -> tuple[torch.Tensor, torch.Tensor]:
    # Take top residue annotations
    top_residue_annotations_idx = logits.argsort(dim=-1, descending=True)[
        ..., :MAX_RESIDUE_ANNOTATIONS
    ]  # (L, MAX_R)
    top_residue_annotations_logprobs = torch.gather(
        F.logsigmoid(logits), -1, top_residue_annotations_idx
    )  # (L, MAX_R)
    top_residue_annotations_probs = top_residue_annotations_logprobs.exp()
    # Keep only positive predictions
    is_negative = top_residue_annotations_probs < annotation_threshold
    top_residue_annotations_idx[is_negative] = 0

    top_residue_annotations_logprobs = top_residue_annotations_logprobs

    return top_residue_annotations_idx, top_residue_annotations_logprobs


def top_p_logits(
    logits: torch.Tensor,
    top_p: float | torch.Tensor,
) -> torch.Tensor:
    top_p = _tensorize_like(top_p, logits)

    batch_dims = logits.size()[:-1]
    logits = logits.reshape(-1, logits.shape[-1])

    # Sort logits in descending order and extract the mask for the top_p
    sorted_logits, sorted_indices = torch.sort(logits, dim=-1, descending=True)
    cumsum_logits = sorted_logits.softmax(-1).cumsum(-1)
    top_p_mask = cumsum_logits <= top_p[:, None]

    # Make sure at least one token is sampled
    top_p_mask[:, 0] = True

    # Mask out the logits that are not in the top_p
    batch_indices_to_mask, _ = torch.where(~top_p_mask)
    vocab_indices_to_mask = sorted_indices[~top_p_mask]
    logits[batch_indices_to_mask, vocab_indices_to_mask] = torch.finfo(logits.dtype).min

    return logits.reshape(*batch_dims, -1)


def _tensorize_like(value: int | float | torch.Tensor, logits: torch.Tensor):
    if isinstance(value, (float, int)):
        value = torch.full_like(logits[..., 0], value, dtype=logits.dtype)
    return value.to(logits.device).expand_as(logits[..., 0]).reshape(-1)
