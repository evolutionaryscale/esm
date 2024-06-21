from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]
PathOrBuffer = Union[PathLike, io.StringIO]


@dataclass
class FunctionAnnotation:
    """Represents an annotation of a protein's function over a range of residues.

    Fields:
        label (str): An entry in either the function_tokens or residue_annotations tokenizer vocabs
        start (int): Start index of this annotation.  1-indexed, inclusive.
        end (int): End index of this annotation.  1-indexed, inclusive.
    """

    label: str
    start: int
    end: int

    def to_tuple(self) -> tuple[str, int, int]:
        return self.label, self.start, self.end
