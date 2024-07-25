from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

PathLike = Union[str, Path, CloudPath]
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

    def __len__(self) -> int:
        """Length of the annotation."""
        return self.end - self.start + 1
