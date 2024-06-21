from esm.tokenization.tokenizer_base import EsmTokenizerBase


class StructureTokenizer(EsmTokenizerBase):
    """A convenince class for accessing special token ids of
    the StructureTokenEncoder and StructureTokenDecoder."""

    def __init__(self, vq_vae_special_tokens: dict[str, int]):
        self.vq_vae_special_tokens = vq_vae_special_tokens

    def mask_token(self) -> str:
        raise NotImplementedError(
            "Structure tokens are defined on 3D coordinates, not strings."
        )

    @property
    def mask_token_id(self) -> int:
        return self.vq_vae_special_tokens["MASK"]

    def bos_token(self) -> str:
        raise NotImplementedError(
            "Structure tokens are defined on 3D coordinates, not strings."
        )

    @property
    def bos_token_id(self) -> int:
        return self.vq_vae_special_tokens["BOS"]

    def eos_token(self) -> str:
        raise NotImplementedError(
            "Structure tokens are defined on 3D coordinates, not strings."
        )

    @property
    def eos_token_id(self) -> int:
        return self.vq_vae_special_tokens["EOS"]

    def pad_token(self) -> str:
        raise NotImplementedError(
            "Structure tokens are defined on 3D coordinates, not strings."
        )

    @property
    def pad_token_id(self) -> int:
        return self.vq_vae_special_tokens["PAD"]

    @property
    def chainbreak_token_id(self) -> int:
        return self.vq_vae_special_tokens["CHAINBREAK"]

    def encode(self, *args, **kwargs):
        raise NotImplementedError(
            "The StructureTokenizer class is provided as a convenience for "
            "accessing special token ids of the StructureTokenEncoder and StructureTokenDecoder.\n"
            "Please use them instead."
        )

    def decode(self, *args, **kwargs):
        raise NotImplementedError(
            "The StructureTokenizer class is provided as a convenience for "
            "accessing special token ids of the StructureTokenEncoder and StructureTokenDecoder.\n"
            "Please use them instead."
        )
