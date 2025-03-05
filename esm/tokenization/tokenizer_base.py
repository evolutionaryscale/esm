from typing import Protocol, runtime_checkable


@runtime_checkable
class EsmTokenizerBase(Protocol):
    mask_token: str
    mask_token_id: int
    bos_token: str
    bos_token_id: int
    eos_token: str
    eos_token_id: int
    pad_token: str
    pad_token_id: int
    chain_break_token: str
    chain_break_token_id: int

    def encode(self, *args, **kwargs): ...

    def decode(self, *args, **kwargs): ...

    @property
    def all_token_ids(self): ...

    @property
    def special_token_ids(self): ...
