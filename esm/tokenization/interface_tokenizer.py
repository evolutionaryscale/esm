from functools import cached_property

import torch

from esm.tokenization.tokenizer_base import EsmTokenizerBase


class InterfaceAnnotationsTokenizer(EsmTokenizerBase):
    """Tokenizer for Interface annotations."""

    @cached_property
    def special_tokens(self) -> list[str]:
        return ["<pad>", "<mask>"]

    @cached_property
    def vocab(self):
        """Tokenizer vocabulary list."""
        nonspecial_tokens = ["<non-interface>", "<interface>"]
        return [*self.special_tokens, *nonspecial_tokens]

    @cached_property
    def vocab_to_index(self) -> dict[str, int]:
        """Constructs token -> token id mapping."""
        return {word: i for i, word in enumerate(self.vocab)}

    def get_special_tokens_mask(self, token_ids: torch.Tensor) -> torch.Tensor:
        """Determines which positions are special tokens.

        Args:
            token_ids: <int>[length]
        Returns:
            <bool>[length] tensor, true where special tokens are located in the input.
        """
        return token_ids < len(self.special_tokens)

    def encode(
        self, annotations: list[str], add_special_tokens: bool = True
    ) -> torch.Tensor:
        """Encode interface annotation list

        Args:
            list[str]: interface annotation e.g. ["<pad>", "<interface>", "<non-interface>"]
        Returns:
            <int>[sequence_length] token ids representing. Will add <bos>/<eos>.
        """
        ids = []
        if add_special_tokens:
            ids.append(self.vocab_to_index["<pad>"])  # bos
        for annot in annotations:
            ids.append(self.vocab_to_index[annot])
        if add_special_tokens:
            ids.append(self.vocab_to_index["<pad>"])  # eos
        return torch.tensor(ids, dtype=torch.int64)

    def decode(self, encoded: torch.Tensor) -> list[str]:
        """Decodes token ids into list of interface annotations.

        Args:
            encoded: <int>[length] token id array.
        Returns
            Decoded secondary structure string.
        """
        return [self.vocab[i] for i in encoded]

    @property
    def mask_token(self) -> str:
        return "<mask>"

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
