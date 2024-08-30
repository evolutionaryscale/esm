from typing import Protocol, runtime_checkable


@runtime_checkable
class EsmTokenizerBase(Protocol):
    def encode(self, *args, **kwargs):
        ...

    def decode(self, *args, **kwargs):
        ...

    @property
    def mask_token(self) -> str:
        ...

    @property
    def mask_token_id(self) -> int:
        ...

    @property
    def bos_token(self) -> str:
        ...

    @property
    def bos_token_id(self) -> int:
        ...

    @property
    def eos_token(self) -> str:
        ...

    @property
    def eos_token_id(self) -> int:
        ...

    @property
    def pad_token(self) -> str:
        ...

    @property
    def pad_token_id(self) -> int:
        ...

    @property
    def chain_break_token(self) -> str:
        ...

    @property
    def chain_break_token_id(self) -> int:
        ...

    @property
    def all_token_ids(self):
        ...

    @property
    def special_token_ids(self):
        ...
