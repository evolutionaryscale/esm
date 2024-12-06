from __future__ import annotations

import contextlib

import attr
import torch
import torch.nn as nn
from attr import dataclass

from esm.layers.regression_head import RegressionHead
from esm.layers.transformer_stack import TransformerStack
from esm.sdk.api import (
    ESMCInferenceClient,
    ESMProtein,
    ESMProteinTensor,
    ForwardTrackData,
    LogitsConfig,
    LogitsOutput,
)
from esm.tokenization import EsmSequenceTokenizer
from esm.utils import encoding
from esm.utils.constants.models import ESMC_600M
from esm.utils.decoding import decode_sequence
from esm.utils.misc import stack_variable_length_tensors
from esm.utils.sampling import _BatchedESMProteinTensor


@dataclass
class ESMCOutput:
    sequence_logits: torch.Tensor
    embeddings: torch.Tensor | None
    hidden_states: torch.Tensor | None


class ESMC(nn.Module, ESMCInferenceClient):
    """
    ESMC model implementation.

    Args:
        d_model (int): The dimensionality of the input and output feature vectors.
        n_heads (int): The number of attention heads in the transformer layers.
        n_layers (int): The number of transformer layers.
    """

    def __init__(
        self, d_model: int, n_heads: int, n_layers: int, tokenizer: EsmSequenceTokenizer
    ):
        super().__init__()
        self.embed = nn.Embedding(64, d_model)
        self.transformer = TransformerStack(
            d_model, n_heads, None, n_layers, n_layers_geom=0
        )
        self.sequence_head = RegressionHead(d_model, 64)
        self.tokenizer = tokenizer

    @classmethod
    def from_pretrained(
        cls, model_name: str = ESMC_600M, device: torch.device | None = None
    ) -> ESMC:
        from esm.pretrained import load_local_model

        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = load_local_model(model_name, device=device)
        if device.type != "cpu":
            model = model.to(torch.bfloat16)
        assert isinstance(model, ESMC)
        return model

    @property
    def device(self):
        return next(self.parameters()).device

    @property
    def raw_model(self):
        return self

    def _tokenize(self, sequence: list[str]) -> torch.Tensor:
        pad = self.tokenizer.pad_token_id
        assert pad is not None
        return stack_variable_length_tensors(
            [
                encoding.tokenize_sequence(x, self.tokenizer, add_special_tokens=True)
                for x in sequence
            ],
            constant_value=pad,
        ).to(next(self.parameters()).device)

    def _detokenize(self, sequence: torch.Tensor) -> list[str]:
        pad = self.tokenizer.pad_token_id
        assert pad is not None
        assert sequence.ndim == 2
        return [decode_sequence(x[x != pad][1:-1], self.tokenizer) for x in sequence]

    def forward(
        self,
        sequence_tokens: torch.Tensor | None = None,
        sequence_id: torch.Tensor | None = None,
    ) -> ESMCOutput:
        """
        Performs forward pass through the ESMC model. Check utils to see how to tokenize inputs from raw data.

        Args:
            sequence_tokens (torch.Tensor, optional): The amino acid tokens.
            sequence_id (torch.Tensor, optional): The sequence ID.

        Returns:
            ESMCOutput: The output of the ESMC model.

        """
        if sequence_id is None:
            sequence_id = sequence_tokens == self.tokenizer.pad_token_id

        x = self.embed(sequence_tokens)
        x, _, hiddens = self.transformer(x, sequence_id=sequence_id)
        sequence_logits = self.sequence_head(x)
        output = ESMCOutput(
            sequence_logits=sequence_logits, embeddings=x, hidden_states=hiddens
        )
        return output

    def encode(self, input: ESMProtein) -> ESMProteinTensor:
        input = attr.evolve(input)  # Make a copy
        sequence_tokens = None

        if input.sequence is not None:
            sequence_tokens = self._tokenize([input.sequence])[0]
        return ESMProteinTensor(sequence=sequence_tokens).to(
            next(self.parameters()).device
        )

    def decode(self, input: ESMProteinTensor) -> ESMProtein:
        input = attr.evolve(input)  # Make a copy

        assert input.sequence is not None
        sequence = self._detokenize(input.sequence)[0]

        return ESMProtein(sequence=sequence)

    def logits(
        self,
        input: ESMProteinTensor | _BatchedESMProteinTensor,
        config: LogitsConfig = LogitsConfig(),
    ) -> LogitsOutput:
        if not isinstance(input, _BatchedESMProteinTensor):
            # Create batch dimension if necessary.
            input = _BatchedESMProteinTensor.from_protein_tensor(input)

        device = torch.device(input.device)

        with (
            torch.no_grad(),
            torch.autocast(enabled=True, device_type=device.type, dtype=torch.bfloat16)  # type: ignore
            if device.type == "cuda"
            else contextlib.nullcontext(),
        ):
            output = self.forward(sequence_tokens=input.sequence)

        return LogitsOutput(
            logits=ForwardTrackData(
                sequence=output.sequence_logits if config.sequence else None
            ),
            embeddings=output.embeddings if config.return_embeddings else None,
        )
