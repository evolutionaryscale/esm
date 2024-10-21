from __future__ import annotations

from abc import ABC
from typing import Sequence

import attr
import torch
from attr import asdict, define

import esm.utils.constants.api as C
from esm.tokenization import (
    TokenizerCollectionProtocol,
    get_model_tokenizers,
)
from esm.utils import encoding
from esm.utils.constants.models import ESM3_OPEN_SMALL
from esm.utils.misc import (
    get_chainbreak_boundaries_from_sequence,
)
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.structure.protein_complex import ProteinComplex
from esm.utils.types import (
    FunctionAnnotation,
    PathOrBuffer,
)


class ProteinType(ABC):
    ...


## Basic Types
@define
class ESMProtein(ProteinType):
    # Tracks
    sequence: str | None = None
    secondary_structure: str | None = None
    sasa: list[float | None] | None = None
    function_annotations: list[FunctionAnnotation] | None = None
    coordinates: torch.Tensor | None = None

    # Metrics
    plddt: torch.Tensor | None = None
    ptm: torch.Tensor | None = None


    # When calling EvolutionaryScale API, use this flag to disclose any
    # sequences that may potentially have concerns.
    # Such sequences may not go through standard safety filter for approved users.
    # Reach out if interested in using this.
    potential_sequence_of_concern: bool = False

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

    @classmethod
    def from_protein_complex(
        cls, protein_complex: ProteinComplex, with_annotations: bool = False
    ) -> ESMProtein:
        if with_annotations:
            raise NotImplementedError(
                "Annotations are not supported for ProteinComplex yet."
            )

        return ESMProtein(
            sequence=protein_complex.sequence,
            secondary_structure=None,
            sasa=None,
            function_annotations=None,
            coordinates=torch.tensor(protein_complex.atom37_positions),
        )

    def to_pdb(self, pdb_path: PathOrBuffer) -> None:
        # Note: Will work for single chains as well and produce same pdb file
        protein_complex = self.to_protein_complex().infer_oxygen()
        protein_complex.to_pdb(pdb_path)

    def to_pdb_string(self) -> str:
        protein_chain = self.to_protein_chain()
        return protein_chain.to_pdb_string()

    def to_protein_chain(self) -> ProteinChain:
        if self.coordinates is None:
            raise ValueError("Coordinates are required to convert to a ProteinChain.")
        protein_chain = ProteinChain.from_atom37(
            atom37_positions=self.coordinates.to("cpu").numpy(),
            id=None,
            sequence=None if self.sequence is None else self.sequence.replace("_", "X"),
            chain_id=None,
            entity_id=None,
            residue_index=None,
            insertion_code=None,
            confidence=None
            if self.plddt is None
            else self.plddt.detach().cpu().numpy(),
        )
        return protein_chain

    def to_protein_complex(
        self, copy_annotations_from_ground_truth: ProteinComplex | None = None
    ) -> ProteinComplex:
        assert (
            self.sequence is not None
        ), "ESMProtein must have a sequence to convert to ProteinComplex"
        assert (
            self.coordinates is not None
        ), "ESMProtein must have coordinates to convert to ProteinComplex"
        coords = self.coordinates.to("cpu").numpy()

        chain_boundaries = get_chainbreak_boundaries_from_sequence(self.sequence)
        if copy_annotations_from_ground_truth is not None:
            gt_chains = list(copy_annotations_from_ground_truth.chain_iter())
        else:
            gt_chains = None
        pred_chains = []
        for i, (start, end) in enumerate(chain_boundaries):
            pred_chain = ProteinChain.from_atom37(
                atom37_positions=coords[start:end],
                sequence=self.sequence[start:end],
                chain_id=gt_chains[i].chain_id if gt_chains is not None else None,
                entity_id=gt_chains[i].entity_id if gt_chains is not None else None,
            )
            pred_chains.append(pred_chain)
        return ProteinComplex.from_chains(pred_chains)


@define
class ESMProteinTensor(ProteinType):
    sequence: torch.Tensor | None = None
    structure: torch.Tensor | None = None
    secondary_structure: torch.Tensor | None = None
    sasa: torch.Tensor | None = None
    function: torch.Tensor | None = None
    residue_annotations: torch.Tensor | None = None
    coordinates: torch.Tensor | None = None

    # When calling EvolutionaryScale API, use this flag to disclose any
    # sequences that may potentially have concerns.
    # Such sequences may not go through standard safety filter for approved users.
    # Reach out if interested in using this.
    potential_sequence_of_concern: bool = False

    def _detect_attribute(self, func, msg):
        mapped = {
            k: func(k, v)
            for k, v in asdict(self).items()
            if isinstance(v, torch.Tensor)
        }
        s = set(mapped.values())
        if len(s) <= 0:
            return None
        if len(s) != 1:
            raise ValueError(f"Either no tracks or inconsistent {msg}: {mapped}")
        return next(iter(s))

    def __len__(self) -> int:
        l = self._detect_attribute(lambda _, x: x.size(0), "length")
        return l if l is not None else 0

    @property
    def device(self) -> str | torch.device:
        d = self._detect_attribute(lambda _, x: x.device, "device")
        assert d is not None
        return d

    def to(self, device_or_dtype: str | torch.device | torch.dtype) -> ESMProteinTensor:
        def _to(name):
            v = getattr(self, name)
            if v is not None and isinstance(v, torch.Tensor):
                setattr(self, name, v.to(device_or_dtype))

        for n in attr.fields(ESMProteinTensor):
            _to(n.name)

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


@define
class ESMProteinError(Exception, ProteinType):
    error_code: int  # Error code follows HTTP convention, i.e., 404 NotFoundError, 500 InternalError.
    error_msg: str


## High Level Endpoint Types
@define
class GenerationConfig:
    track: str = ""
    invalid_ids: Sequence[int] = []
    schedule: str = "cosine"
    # Set this to a higher value for better generation results.
    # Note that this needs to be less than or equal to the sequence length.
    num_steps: int = 1
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
    sequence: SamplingTrackConfig | None = attr.field(
        default=None, metadata={"max_topk": C.MAX_TOPK_SEQUENCE}
    )
    structure: SamplingTrackConfig | None = attr.field(
        default=None, metadata={"max_topk": C.MAX_TOPK_STRUCTURE}
    )
    secondary_structure: SamplingTrackConfig | None = attr.field(
        default=None, metadata={"max_topk": C.MAX_TOPK_SECONDARY_STRUCTURE}
    )
    sasa: SamplingTrackConfig | None = attr.field(
        default=None, metadata={"max_topk": C.MAX_TOPK_SASA}
    )
    function: SamplingTrackConfig | None = attr.field(
        default=None, metadata={"max_topk": C.MAX_TOPK_FUNCTION}
    )

    return_per_residue_embeddings: bool = False
    return_mean_embedding: bool = False


@define
class ForwardTrackData:
    sequence: torch.Tensor | None = None
    structure: torch.Tensor | None = None
    secondary_structure: torch.Tensor | None = None
    sasa: torch.Tensor | None = None
    function: torch.Tensor | None = None


@define
class LogitsConfig:
    # Logits.
    sequence: bool = False
    structure: bool = False
    secondary_structure: bool = False
    sasa: bool = False
    function: bool = False
    residue_annotations: bool = False

    # Embeddings.
    return_embeddings: bool = False


@define
class LogitsOutput:
    logits: ForwardTrackData | None = None
    embeddings: torch.Tensor | None = None

    # Residue annotations is multi-hot, so deserves special treatment
    # It's not a categorical distribution, but instead a bernoulli, so
    # softmax across the last dimension is _wrong_
    residue_annotation_logits: torch.Tensor | None = None


@define
class ForwardAndSampleOutput(LogitsOutput):
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


class ESM3InferenceClient(ABC):
    def generate(self, input: ProteinType, config: GenerationConfig) -> ProteinType:
        # This is the easiest and most flexible way to run ESM3. Generate will
        # iteratively sample tokens an provide an output with the track specified
        # completely filled out, according to the GenerationConfig provided.
        # It is a local function wrapping calls for encode -> iterative_sampling -> decode.
        # if a ESMProteinTensor is provided, encode and decode are skipped
        raise NotImplementedError

    def batch_generate(
        self,
        inputs: Sequence[ProteinType],
        configs: Sequence[GenerationConfig],
    ) -> Sequence[ProteinType]:
        # Same as generate(...), but generates a batch of proteins at once.
        raise NotImplementedError

    def encode(self, input: ESMProtein) -> ESMProteinTensor:
        # Encode allows for encoding RawRepresentation into TokenizedRepresentation.
        # This runs the structure_token_encoder, as well as dealing with PDB => atom37 conversion
        raise NotImplementedError

    def decode(self, input: ESMProteinTensor) -> ESMProtein:
        # Decode is the inverse of encode, and runs a structure_token_decoder to output coordinates
        raise NotImplementedError

    def logits(
        self, input: ESMProteinTensor, config: LogitsConfig = LogitsConfig()
    ) -> LogitsOutput:
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

    @property
    def raw_model(self):
        # Get underlying esm3 model of an inference client.
        raise NotImplementedError
