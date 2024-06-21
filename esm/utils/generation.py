from typing import Callable

import attr
import torch
from tqdm import tqdm

from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    ESMProteinTensor,
    GenerationConfig,
    SamplingConfig,
    SamplingTrackConfig,
)
from esm.tokenization import (
    EsmTokenizerBase,
    TokenizerCollectionProtocol,
)
from esm.utils.constants import esm3 as C
from esm.utils.noise_schedules import NOISE_SCHEDULE_REGISTRY


def iterative_sampling_raw(
    client: ESM3InferenceClient,
    input: ESMProtein,
    config: GenerationConfig,
):
    # Keep structure tokens
    input_tokens = client.encode(input)

    output_tokens = client.generate(input_tokens, config)

    raw_protein = client.decode(output_tokens)

    track_to_sample = config.track

    if track_to_sample not in ["function", "residue_annotations"]:
        # Function and residue annotation encoding/decoding is lossy
        # There is no guarantee that decoding encoded tokens will yield the same input
        raw_protein.function_annotations = input.function_annotations

    return raw_protein


def iterative_sampling_tokens(
    client: ESM3InferenceClient,
    input_tokens: ESMProteinTensor,
    config: GenerationConfig,
    tokenizers: TokenizerCollectionProtocol,
) -> ESMProteinTensor:
    track_to_sample = config.track

    # Get all tracks that require sampling
    all_tracks = [
        f.name for f in attr.fields(SamplingConfig) if "embedding" not in f.name
    ]

    sequence_length = len(input_tokens)
    device = input_tokens.device

    # Initialize schedule and masks
    decoding_schedule = NOISE_SCHEDULE_REGISTRY[config.schedule]
    sampled_tokens = attr.evolve(input_tokens)  # Make a copy

    if config.condition_on_coordinates_only and input_tokens.coordinates is not None:
        sampled_tokens.structure = None

    sampling_mask = torch.ones(
        sequence_length,
        dtype=torch.bool,
        device=device,
    )
    sampling_mask[0] = False
    sampling_mask[-1] = False

    get_tokenizer: Callable[[str], EsmTokenizerBase] = lambda s: getattr(tokenizers, s)
    if getattr(sampled_tokens, track_to_sample) is None:
        if track_to_sample == "function":
            dims = (sequence_length, tokenizers.function.depth)
        elif track_to_sample == "residue_annotations":
            dims = (sequence_length, C.MAX_RESIDUE_ANNOTATIONS)
        else:
            dims = (sequence_length,)
        masked_tokens = torch.full(
            dims,
            get_tokenizer(track_to_sample).mask_token_id,
            dtype=torch.long,
            device=device,
        )
        if track_to_sample == "sequence":
            masked_tokens[0] = tokenizers.sequence.cls_token_id  # type: ignore
            masked_tokens[-1] = tokenizers.sequence.eos_token_id  # type: ignore
        else:
            masked_tokens[0] = get_tokenizer(track_to_sample).bos_token_id
            masked_tokens[-1] = get_tokenizer(track_to_sample).eos_token_id

        setattr(
            sampled_tokens,
            track_to_sample,
            masked_tokens,
        )
    else:
        is_mask: torch.Tensor = (
            getattr(input_tokens, track_to_sample)
            == get_tokenizer(track_to_sample).mask_token_id
        )
        if not is_mask.any().item():
            raise ValueError(f"Cannot sample {config.track} when input has no masks.")
        sampling_mask = sampling_mask & is_mask

    # Decode

    def maybe_clone(x: torch.Tensor | None) -> torch.Tensor | None:
        return x.clone() if x is not None else None

    L = sequence_length - 2
    positions_sampled = 0
    for t in tqdm(range(config.num_steps)):
        # Single step sampling at all positions
        track_sample_config = SamplingTrackConfig()
        track_sample_config.invalid_ids = config.invalid_ids
        track_sample_config.temperature = config.temperature
        track_sample_config.top_p = config.top_p
        sampling_config = SamplingConfig(**{track_to_sample: track_sample_config})  # type: ignore

        forward_and_sample_output = client.forward_and_sample(
            sampled_tokens, sampling_config
        )
        new_samples = forward_and_sample_output.protein_tensor

        # Calculate number of tokens to sample
        perc_masked = decoding_schedule(torch.tensor((t + 1) / config.num_steps))
        num_to_sample = int((1 - perc_masked) * L) - positions_sampled
        positions_sampled += num_to_sample

        # Select tokens based on lowest entropy
        if track_to_sample in ["function", "residue_annotations"]:
            # TODO: Implement iterative decoding for function and residue_annotations
            # TODO: Fix encode/decode of interpro tokens (not yet supported)
            sampled_tokens.function = maybe_clone(input_tokens.function)
            sampled_tokens.residue_annotations = maybe_clone(
                input_tokens.residue_annotations
            )
            if track_to_sample in track_to_sample:
                raise NotImplementedError(
                    f"Iterative decoding for {track_to_sample} is not supported yet."
                )
            continue

        sampling_mask = sampling_mask & (
            getattr(sampled_tokens, track_to_sample)
            == get_tokenizer(track_to_sample).mask_token_id
        )

        track_entropy: torch.Tensor = getattr(
            forward_and_sample_output.entropy, track_to_sample
        )
        track_entropy = track_entropy.masked_fill(
            ~sampling_mask, torch.finfo(track_entropy.dtype).max
        )
        _, indices = track_entropy.topk(num_to_sample, dim=-1, largest=False)
        is_top_k = ~(
            torch.arange(sequence_length, device=device)[:, None] != indices[None, :]
        ).all(-1)
        tokens_to_sample = sampling_mask & is_top_k

        old_track_samples = getattr(sampled_tokens, track_to_sample)
        new_track_samples = getattr(new_samples, track_to_sample)

        new_track_samples = torch.where(
            tokens_to_sample, new_track_samples, old_track_samples
        )

        setattr(sampled_tokens, track_to_sample, new_track_samples)

    # Do not update tracks that were not sampled (e.g. keep None instead of masks)
    for track in all_tracks:
        if track != track_to_sample:
            setattr(
                sampled_tokens,
                track,
                maybe_clone(getattr(input_tokens, track)),
            )

    return sampled_tokens
