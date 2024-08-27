from functools import cached_property
from typing import Sequence

import torch

from esm.tokenization.tokenizer_base import EsmTokenizerBase
from esm.utils.constants import esm3 as C


class SecondaryStructureTokenizer(EsmTokenizerBase):
    """Tokenizer for secondary structure strings."""

    def __init__(self, kind: str = "ss8"):
        assert kind in ("ss8", "ss3")
        self.kind = kind

    @property
    def special_tokens(self) -> list[str]:
        return ["<pad>", "<motif>", "<unk>"]

    @cached_property
    def vocab(self):
        """Tokenzier vocabulary list."""
        match self.kind:
            case "ss8":
                nonspecial_tokens = list(C.SSE_8CLASS_VOCAB)  # "GHITEBSC"
            case "ss3":
                nonspecial_tokens = list(C.SSE_3CLASS_VOCAB)  # HEC
            case _:
                raise ValueError(self.kind)

        # The non-special tokens ids match amino acid tokens ids when possible.
        return [*self.special_tokens, *nonspecial_tokens]

    @cached_property
    def vocab_to_index(self) -> dict[str, int]:
        """Constructs token -> token id mapping."""
        return {word: i for i, word in enumerate(self.vocab)}

    def get_special_tokens_mask(self, tokens: torch.Tensor) -> torch.Tensor:
        """Determines which positions are special tokens.

        Args:
            tokens: <int>[length]
        Returns:
            <bool>[length] tensor, true where special tokens are located in the input.
        """
        return tokens < len(self.special_tokens)

    def encode(
        self, sequence: str | Sequence[str], add_special_tokens: bool = True
    ) -> torch.Tensor:
        """Encode secondary structure string

        Args:
            string: secondary structure string e.g. "GHHIT", or as token listk.
        Returns:
            <int>[sequence_length] token ids representing. Will add <cls>/<eos>.
        """
        ids = []
        if add_special_tokens:
            ids.append(self.vocab_to_index["<pad>"])  # cls
        for char in sequence:
            ids.append(self.vocab_to_index[char])
        if add_special_tokens:
            ids.append(self.vocab_to_index["<pad>"])  # eos
        return torch.tensor(ids, dtype=torch.int64)

    def decode(self, encoded: torch.Tensor) -> str:
        """Decodes token ids into secondary structure string.

        Args:
            encoded: <int>[length] token id array.
        Returns
            Decoded secondary structure string.
        """
        return "".join(self.vocab[i] for i in encoded)

    @property
    def mask_token(self) -> str:
        return "<pad>"

    @property
    def mask_token_id(self) -> int:
        return self.vocab_to_index[self.mask_token]

    @property
    def bos_token(self) -> str:
        return "<pad>"

    @property
    def bos_token_id(self) -> int:
        return self.vocab_to_index[self.bos_token]

    @property
    def eos_token(self) -> str:
        return "<pad>"

    @property
    def eos_token_id(self) -> int:
        return self.vocab_to_index[self.eos_token]

    @property
    def pad_token(self) -> str:
        return "<pad>"

    @property
    def pad_token_id(self) -> int:
        return self.vocab_to_index[self.pad_token]

    @property
    def chain_break_token(self) -> str:
        return "<pad>"

    @property
    def chain_break_token_id(self) -> int:
        return self.vocab_to_index[self.chain_break_token]

    @property
    def all_token_ids(self):
        return list(range(len(self.vocab)))

    @property
    def special_token_ids(self):
        return [self.vocab_to_index[token] for token in self.special_tokens]
