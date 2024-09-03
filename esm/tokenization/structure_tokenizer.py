from esm.tokenization.tokenizer_base import EsmTokenizerBase
from esm.utils.constants import esm3 as C


class StructureTokenizer(EsmTokenizerBase):
    """A convenince class for accessing special token ids of
    the StructureTokenEncoder and StructureTokenDecoder."""

    def __init__(self, codebook_size: int = C.VQVAE_CODEBOOK_SIZE):
        self.vq_vae_special_tokens = {
            "MASK": codebook_size,
            "EOS": codebook_size + 1,
            "BOS": codebook_size + 2,
            "PAD": codebook_size + 3,
            "CHAINBREAK": codebook_size + 4,
        }

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

    def chain_break_token(self) -> str:
        raise NotImplementedError(
            "Structure tokens are defined on 3D coordinates, not strings."
        )

    @property
    def chain_break_token_id(self) -> int:
        return self.vq_vae_special_tokens["CHAINBREAK"]

    @property
    def all_token_ids(self):
        return list(range(C.VQVAE_CODEBOOK_SIZE + len(self.vq_vae_special_tokens)))

    @property
    def special_token_ids(self):
        return self.vq_vae_special_tokens.values()

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
