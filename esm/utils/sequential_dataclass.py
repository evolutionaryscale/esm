from abc import ABC, abstractmethod
from dataclasses import dataclass, fields, replace
from typing import TypeVar

import numpy as np

from esm.utils.misc import concat_objects, slice_any_object

T = TypeVar("T")


@dataclass(frozen=True)
class SequentialDataclass(ABC):
    """
    This is a builder on a dataclass that allows for automatic slicing and concatenation.

    When representing multimodal data, we often have multiple datatypes which have sequence dimensions that are the same (e.g. the length of the protein).

    When applying a transformation like a crop, we want to apply this to all tensors at the same time (e.g. crop the sequence, structure, and function).

    We also have some fields that are not sequential (like an id, or data source), which we don't want to crop.

    The SequentialDataclass abstracts this cropping away, allowing you to define dataclasses that implement `__len__`, `__getitem__` and `concat` automatically.

    This is done through the `metadata` field, which can take 3 values:
        `sequence` (bool): True or False, tells the dataclass whether this field is a sequential type. Default: False.
        `sequence_dim` (int): Which dimension is the sequential dimension (e.g. for a list of inverse folded sequences, we want to index each sequence in the list, not the list itself). Default: 0.
        `join_token` (Any): What token to use to join when concatenating elements. Default: None.


    Example:

        @dataclass(frozen=True)
        class Foo(SequentialDataclass):
            id: str
            sequence: str = field(metadata={"sequence": True, "join_token": "|"})
            tensor: torch.Tensor = field(metadata={"sequence": True, "join_token": torch.nan})

            def __len__(self):
                # Must implement the __len__ method
                return len(self.sequence)

        >>> foo = Foo(id="foo", sequence="ABCDE", tensor=torch.randn(5))
        Foo(id='foo', sequence='ABCDE', tensor=tensor([ 0.0252, -0.3335, -0.5143,  0.0251, -1.0717]))

        >>> foo[1:4]
        Foo(id='foo', sequence='BCD', tensor=tensor([-0.3335, -0.5143,  0.0251]))

        >>> foo[np.arange(5) < 3]
        Foo(id='foo', sequence='ABC', tensor=tensor([ 0.0252, -0.3335, -0.5143]))

        >>> Foo.concat([foo[:2], foo[3:]])
        Foo(id='foo', sequence='AB|DE', tensor=tensor([ 0.0252, -0.3335,     nan,  0.0251, -1.0717]))

        # Trying to create a type where the sequence lengths do not match raises an error
        >>> foo = Foo(id="foo", sequence="ABCDE", tensor=torch.randn(6))
        ValueError: Mismatch in sequence length for field: tensor. Expected 5, received 6

    """

    def __post_init__(self):
        self._check_sequence_lengths_match()

    @abstractmethod
    def __len__(self):
        raise NotImplementedError

    def __getitem__(self, idx: int | list[int] | slice | np.ndarray):
        updated_fields = {}
        if isinstance(idx, int):
            # make it so that things remain sequential
            idx = [idx]

        for fld in fields(self):
            if fld.metadata.get("sequence", False):
                # this is a sequence, should be the same length as all other sequences
                sequence_dim = fld.metadata.get("sequence_dim", 0)
                value = getattr(self, fld.name)
                if value is None:
                    continue
                match sequence_dim:
                    case 0:
                        # sequence is first dimension
                        value = getattr(self, fld.name)
                        value = slice_any_object(value, idx)
                        updated_fields[fld.name] = value
                    case 1:
                        new_value = [slice_any_object(item, idx) for item in value]
                        updated_fields[fld.name] = value.__class__(new_value)
                    case _:
                        raise NotImplementedError(
                            "Arbitrary slicing for different sequence length fields is not implemented"
                        )

        return replace(self, **updated_fields)

    def _check_sequence_lengths_match(self):
        """Checks if sequence lengths of all "sequence" fields match."""
        for fld in fields(self):
            if fld.metadata.get("sequence", False) and fld.name != "complex":
                # this is a sequence, should be the same length as all other sequences
                sequence_dim = fld.metadata.get("sequence_dim", 0)
                value = getattr(self, fld.name)
                if value is None:
                    continue
                match sequence_dim:
                    case 0:
                        # sequence is first dimension
                        value = getattr(self, fld.name)
                        if len(value) != len(self):
                            raise ValueError(
                                f"Mismatch in sequence length for field: {fld.name}. Expected {len(self)}, received {len(value)}"
                            )
                    case 1:
                        for item in value:
                            if len(item) != len(self):
                                raise ValueError(
                                    f"Mismatch in sequence length for field: {fld.name}. Expected {len(self)}, received {len(item)}"
                                )
                    case _:
                        raise NotImplementedError(
                            "Arbitrary matching for different sequence length fields is not implemented"
                        )

    @classmethod
    def concat(cls, items: list[T], **kwargs) -> T:
        updated_fields = {}
        for fld in fields(cls):
            if fld.metadata.get("sequence", False):
                # this is a sequence, should be the same length as all other sequences
                sequence_dim = fld.metadata.get("sequence_dim", 0)
                join_value = fld.metadata.get("join_token", None)
                if getattr(items[0], fld.name) is None:
                    continue
                values = [getattr(item, fld.name) for item in items]
                match sequence_dim:
                    case 0:
                        # sequence is first dimension
                        value = concat_objects(values, join_value)
                        updated_fields[fld.name] = value
                    case 1:
                        new_value = [
                            concat_objects(item, join_value) for item in zip(*values)
                        ]
                        updated_fields[fld.name] = getattr(
                            items[0], fld.name
                        ).__class__(new_value)
                    case _:
                        raise NotImplementedError(
                            "Arbitrary joining for different sequence length fields is not implemented"
                        )
        updated_fields.update(kwargs)

        return replace(
            items[0],  # type: ignore
            **updated_fields,
        )
