from __future__ import annotations

from abc import ABC
from typing import Sequence, TypeVar

import attr
import torch
from attr import define

from esm.tokenization import (
    TokenizerCollectionProtocol,
    get_model_tokenizers,
)
from esm.utils import encoding
from esm.utils.constants.models import ESM3_OPEN_SMALL
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.types import (
    FunctionAnnotation,
    PathLike,
    PathOrBuffer,
)


## Basic Types
@define
class ESMProtein:
    # Tracks
    sequence: str | None = None
    secondary_structure: str | None = None
    sasa: list[float | str | None] | None = None
    function_annotations: list[FunctionAnnotation] | None = None
    coordinates: torch.Tensor | None = None
    # Metrics
    plddt: torch.Tensor | None = None
    ptm: torch.Tensor | None = None

    def __len__(self):
        if self.sequence is not None:
            return len(self.sequence)
        elif self.secondary_structure is not None:
            return len(self.secondary_structure)
        elif self.sasa is not None:
            return len(self.sasa)
        elif self.coordinates is not None:
            return self.coordinates.size(0)
        else:
            raise ValueError("No track to determine length from.")

    @classmethod
    def from_pdb(
        cls,
        path: PathOrBuffer,
        chain_id: str = "detect",
        id: str | None = None,
        is_predicted: bool = False,
    ) -> ESMProtein:
        protein_chain = ProteinChain.from_pdb(
            path=path, chain_id=chain_id, id=id, is_predicted=is_predicted
        )
        return cls.from_protein_chain(protein_chain)

    @classmethod
    def from_protein_chain(
        cls, protein_chain: ProteinChain, with_annotations: bool = False
    ) -> ESMProtein:
        # By default, we don't annotate with DSSP / SASA, which are expensive.
        # If mkdssp is installed, we can annotate with a flag.
        if with_annotations:
            return ESMProtein(
                sequence=protein_chain.sequence,
                secondary_structure=protein_chain.dssp().tolist(),
                sasa=protein_chain.sasa().tolist(),
                function_annotations=None,
                coordinates=torch.tensor(protein_chain.atom37_positions),
            )
        else:
            return ESMProtein(
                sequence=protein_chain.sequence,
                secondary_structure=None,
                sasa=None,
                function_annotations=None,
                coordinates=torch.tensor(protein_chain.atom37_positions),
            )

    def to_pdb(self, pdb_path: PathLike) -> None:
        protein_chain = self.to_protein_chain()
        protein_chain.to_pdb(pdb_path)

    def to_pdb_string(self) -> str:
        protein_chain = self.to_protein_chain()
        return protein_chain.to_pdb_string()

    def to_protein_chain(self) -> ProteinChain:
        if self.coordinates is None:
            raise ValueError("Coordinates are required to convert to a ProteinChain.")
        protein_chain = ProteinChain.from_atom37(
            atom37_positions=self.coordinates.to("cpu").numpy(),
            id=None,
            sequence=self.sequence,
            chain_id=None,
            entity_id=None,
            residue_index=None,
            insertion_code=None,
            confidence=None if self.plddt is None else self.plddt.detach().cpu().numpy(),
        )
        return protein_chain


@define
class ESMProteinTensor:
    sequence: torch.Tensor | None = None
    structure: torch.Tensor | None = None
    secondary_structure: torch.Tensor | None = None
    sasa: torch.Tensor | None = None
    function: torch.Tensor | None = None
    residue_annotations: torch.Tensor | None = None
    coordinates: torch.Tensor | None = None

    def __len__(self) -> int:
        if self.sequence is not None:
            return self.sequence.size(0)
        elif self.structure is not None:
            return self.structure.size(0)
        elif self.secondary_structure is not None:
            return self.secondary_structure.size(0)
        elif self.sasa is not None:
            return self.sasa.size(0)
        elif self.coordinates is not None:
            return self.coordinates.size(0)
        else:
            raise ValueError("No track to determine length from.")

    @property
    def device(self) -> str | torch.device:
        device_ = None

        tracks = [f.name for f in attr.fields(ESMProteinTensor)]

        for track in tracks:
            current_track: torch.Tensor | None = getattr(self, track)
            if current_track is not None:
                if device_ is not None and device_ != current_track.device:
                    raise ValueError(f"Inconsistent devices for track {track}.")
                device_ = getattr(self, track).device

        if device_ is None:
            raise ValueError("No track to determine device from.")

        return device_

    def to(self, device: str | torch.device | None) -> ESMProteinTensor:
        if device is None:
            return self

        device = torch.device(device)

        def _to(name):
            v = getattr(self, name)
            if v is not None:
                setattr(self, name, v.to(device))

        for n in [
            "sequence",
            "structure",
            "secondary_structure",
            "sasa",
            "function",
            "residue_annotations",
            "coordinates",
        ]:
            _to(n)

        return self

    @classmethod
    def empty(
        cls,
        length: int,
        tokenizers: TokenizerCollectionProtocol | None = None,
        device: torch.device | str = "cpu",
    ) -> ESMProteinTensor:
        if tokenizers is None:
            tokenizers = get_model_tokenizers(ESM3_OPEN_SMALL)

        return ESMProteinTensor(
            sequence=encoding.get_default_sequence_tokens(
                length, tokenizers.sequence
            ).to(device),
            structure=encoding.get_default_structure_tokens(
                length, tokenizers.structure
            ).to(device),
            secondary_structure=encoding.get_default_secondary_structure_tokens(
                length, tokenizers.secondary_structure
            ).to(device),
            sasa=encoding.get_default_sasa_tokens(length, tokenizers.sasa).to(device),
            function=encoding.get_default_function_tokens(
                length, tokenizers.function
            ).to(device),
            residue_annotations=encoding.get_default_residue_annotation_tokens(
                length, tokenizers.residue_annotations
            ).to(device),
        )


## High Level Endpoint Types
@define
class GenerationConfig:
    track: str = ""
    invalid_ids: Sequence[int] = []
    schedule: str = "cosine"
    num_steps: int = 8
    temperature: float = 1.0
    top_p: float = 1.0
    condition_on_coordinates_only: bool = True


## Low Level Endpoint Types
@define
class SamplingTrackConfig:
    temperature: float = 1.0
    top_p: float = 1.0
    only_sample_masked_tokens: bool = True
    invalid_ids: Sequence[int] = []
    topk_logprobs: int = 0


@define
class SamplingConfig:
    sequence: SamplingTrackConfig | None = None
    structure: SamplingTrackConfig | None = None
    secondary_structure: SamplingTrackConfig | None = None
    sasa: SamplingTrackConfig | None = None
    function: SamplingTrackConfig | None = None

    return_per_residue_embeddings: bool = False
    return_mean_embedding: bool = False


@define
class ReturnLogitsConfig:
    sequence: bool = False
    structure: bool = False
    secondary_structure: bool = False
    sasa: bool = False
    function: bool = False
    residue_annotations: bool = False


@define
class ForwardConfig:
    return_logits: ReturnLogitsConfig = ReturnLogitsConfig()
    return_embeddings: bool = False


@define
class ForwardTrackData:
    sequence: torch.Tensor | None = None
    structure: torch.Tensor | None = None
    secondary_structure: torch.Tensor | None = None
    sasa: torch.Tensor | None = None
    function: torch.Tensor | None = None


@define
class ForwardOutput:
    logits: ForwardTrackData | None = None
    embeddings: torch.Tensor | None = None

    # Residue annotations is multi-hot, so deserves special treatment
    # It's not a categorical distribution, but instead a bernoulli, so
    # softmax across the last dimension is _wrong_
    residue_annotation_logits: torch.Tensor | None = None


@define
class ForwardAndSampleOutput(ForwardOutput):
    protein_tensor: ESMProteinTensor = ESMProteinTensor()

    entropy: ForwardTrackData | None = None
    # Probability of sampled token
    prob: ForwardTrackData | None = None
    logprob: ForwardTrackData | None = None
    # Top probability at this position
    top_prob: ForwardTrackData | None = None
    topk_logprob: ForwardTrackData | None = None
    # Which tokens correspond to top probability
    topk_tokens: ForwardTrackData | None = None

    per_residue_embedding: torch.Tensor | None = None
    mean_embedding: torch.Tensor | None = None


ProteinType = TypeVar("ProteinType", bound=ESMProteinTensor | ESMProtein)


class ESM3InferenceClient(ABC):
    def generate(self, input: ProteinType, config: GenerationConfig) -> ProteinType:
        # This is the easiest and most flexible way to run ESM3. Generate will
        # iteratively sample tokens an provide an output with the track specified
        # completely filled out, according to the GenerationConfig provided.
        # It is a local function wrapping calls for encode -> iterative_sampling -> decode.
        # if a ESMProteinTensor is provided, encode and decode are skipped
        raise NotImplementedError

    def encode(self, input: ESMProtein) -> ESMProteinTensor:
        # Encode allows for encoding RawRepresentation into TokenizedRepresentation.
        # This runs the structure_token_encoder, as well as dealing with PDB => atom37 conversion
        raise NotImplementedError

    def decode(self, input: ESMProteinTensor) -> ESMProtein:
        # Decode is the inverse of encode, and runs a structure_token_decoder to output coordinates
        raise NotImplementedError

    def _forward(
        self, input: ESMProteinTensor, config: ForwardConfig = ForwardConfig()
    ) -> ForwardOutput:
        # Our API generally discourages using raw forwards.
        # This is because sending logits can be prohibitively expensive.
        # Please use forward_and_sample instead.
        raise NotImplementedError

    def forward_and_sample(
        self, input: ESMProteinTensor, sampling_configuration: SamplingConfig
    ) -> ForwardAndSampleOutput:
        # forward_and_sample runs a single model forward, sampling tokens according to `SamplingConfiguration`.
        # This is the way for power users to run ESM3. We hope to design this in a way to enable high throughput
        # inference, as well as arbitrary chain-of-though invocations of ESM3.
        raise NotImplementedError
