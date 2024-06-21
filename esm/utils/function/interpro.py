"""Utilities for interacting with InterPro."""

import itertools
import re
from dataclasses import dataclass
from enum import IntEnum, auto
from functools import cached_property
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from esm.utils.constants import esm3 as C


def parse_go_terms(text: str) -> list[str]:
    """Parses GO terms from a string.

    Args:
        text: String containing GO terms. Example: "GO:0008309, GO:1902267" Note that GO
          terms have exactly 7 digits.
    Returns:
        All GO terms found in the string. Example: ['GO:0008309', 'GO:1902267']
    """
    return re.findall(r"GO:(?:\d{7,})", text)


def _parse_interpro2go(path: str) -> dict[str, list[str]]:
    """Parses InterPro2GO file into map.

    NOTE: this file has a very strange, non-standard format.

    Args:
        path: path to InterPro2GO file from: https://www.ebi.ac.uk/GOA/InterPro2GO
    Returns:
        Mapping from InterPro to list of associated GO terms.
    """
    with Path(path).open("r") as f:
        text = f.read()
    df = pd.Series(text.split("\n"), name="line").to_frame()
    df = df[~df.line.str.startswith("!")]
    df["interpro_id"] = df.line.apply(lambda line: re.findall(r"IPR\d+", line))
    df["go_ids"] = df.line.apply(parse_go_terms)
    df = df[df.go_ids.apply(len).gt(0) & df.interpro_id.apply(len).eq(1)]
    df["interpro_id"] = df["interpro_id"].apply(lambda xs: xs[0])  # type: ignore

    # Group all mappints together into a single map.
    df = (
        df.groupby("interpro_id")["go_ids"]  # type: ignore
        .apply(lambda group: list(itertools.chain.from_iterable(group)))
        .reset_index()
    )
    return dict(zip(df.interpro_id, df.go_ids))  # type: ignore


class InterProEntryType(IntEnum):
    """InterPro types and representation counts:

    Family                    21,942
    Domain                    14,053
    Homologous_superfamily     3,446
    Conserved_site               728
    Repeat                       374
    Active_site                  133
    Binding_site                  75
    PTM                           17
    """

    ACTIVE_SITE = 0
    BINDING_SITE = auto()
    CONSERVED_SITE = auto()
    DOMAIN = auto()
    FAMILY = auto()
    HOMOLOGOUS_SUPERFAMILY = auto()
    PTM = auto()
    REPEAT = auto()
    UNKNOWN = auto()


@dataclass
class InterProEntry:
    """Represents an InterPro entry."""

    id: str  # Example: IPR000006
    type: InterProEntryType
    name: str  # Example: "Metallothionein, vertebrate"
    description: str | None = None


@dataclass(frozen=True)
class InterProRangeAnnotation:
    """Represents a InterPro annotation along a range of residues in a protein."""

    interpro_accession: str
    start_idx: int
    end_idx: int


class InterPro:
    """Convenience class interacting with InterPro ontology/data."""

    def __init__(
        self,
        entries_path: str | None = None,
        hierarchy_path: str | None = None,
        interpro2go_path: str | None = None,
    ):
        """Constructs interface to query InterPro entries."""
        default = lambda x, d: x if x is not None else d
        self.entries_path = default(entries_path, str(C.data_root() / C.INTERPRO_ENTRY))
        self.hierarchy_graph_path = default(
            hierarchy_path, str(C.data_root() / C.INTERPRO_HIERARCHY)
        )
        self.interpro2go_path = default(
            interpro2go_path, str(C.data_root() / C.INTERPRO2GO)
        )

    @cached_property
    def interpro2go(self) -> dict[str, list[str]]:
        """Reads the InterPro to GO term mapping."""
        assert self.interpro2go_path is not None
        return _parse_interpro2go(self.interpro2go_path)

    @cached_property
    def entries_frame(self) -> pd.DataFrame:
        """Loads full InterPro entry set as a DataFrame.

        Colums are
            - "id": str interpro accession /id as
            - "type": InterProEntryType representing the type of annotation.
            - "name": Short name of the entry.
        """
        with Path(self.entries_path).open("r") as f:
            df = pd.read_csv(f, sep="\t")
        assert all(
            col in df.columns for col in ["ENTRY_AC", "ENTRY_TYPE", "ENTRY_NAME"]
        )
        df.rename(
            columns={
                "ENTRY_AC": "id",
                "ENTRY_TYPE": "type",
                "ENTRY_NAME": "name",
            },
            inplace=True,
        )
        df["type"] = df.type.str.upper().apply(
            lambda type_name: InterProEntryType[type_name]
        )
        return df

    @cached_property
    def entries(self) -> dict[str, InterProEntry]:
        """Returns all InterPro entries."""
        return {
            row.id: InterProEntry(  # type: ignore
                id=row.id,  # type: ignore
                type=row.type,  # type: ignore
                name=row.name,  # type: ignore
            )
            for row in self.entries_frame.itertuples()
        }

    def lookup_name(self, interpro_id: str) -> str | None:
        """Short name / title for an interpro id."""
        if interpro_id not in self.entries:
            return None
        return self.entries[interpro_id].name

    def lookup_entry_type(self, interpro_id: str) -> InterProEntryType:
        """Looks up entry-type for an interpro id."""
        if interpro_id in self.entries:
            return self.entries[interpro_id].type
        else:
            return InterProEntryType.UNKNOWN

    @cached_property
    def graph(self) -> nx.DiGraph:
        """Reads the InterPro hierarchy of InterPro."""
        graph = nx.DiGraph()
        with Path(self.hierarchy_graph_path).open("r") as f:
            parents = []
            for line in f:
                ipr = line.split("::", maxsplit=1)[0]
                ipr_strip = ipr.lstrip("-")
                level = (len(ipr) - len(ipr_strip)) // 2
                parents = parents[:level]
                graph.add_node(ipr_strip)
                if parents:
                    graph.add_edge(ipr_strip, parents[-1])
                parents.append(ipr_strip)
        return graph


def parse_interpro_features(
    interpro_accessions: list[str],
    interpro_starts: list[int],
    interpro_ends: list[int],
) -> list[InterProRangeAnnotation]:
    """Parses raw InterPro ranges.

    Args:
        interpro_accessions: list of InterPro accessions
        interpro_starts: list of one-indexed inclusive residue locations where the
          annotation from `interpro_accesisons` begin.
        interpro_ends: list of one-indexed *inclusive* residue locations where the
          annotation from `interpro_accesisons` end.
    Returns:
        Collated InterProRangeAnnotations. NOTE that index conversion will convert range
        bounds to zero-indexed [inclusive, exclusive) start/end indices.
    """
    assert len(interpro_accessions) == len(interpro_starts) == len(interpro_ends)

    # Residue locations from Uniprot/InterPro are [inclusive, inclusive] and 1-index.
    start_idcs = np.array(interpro_starts).astype(int)
    end_idcs = np.array(interpro_ends).astype(int)

    # We want to use Python's convention of [inclusive, exclusive) and 0-indexing.
    # Interpro residue indices are [inclusive, inclusive] and 1-indexing.
    # The conversion ends up being:
    #   ```python
    #   end_idcs += 1  # [inclusive, inclusive] -> [inclusive, exclusive)
    #   start_idcs -= 1  # 1 -> 0 indexing
    #   end_idcs -= 1  # 1 -> 0 indexing
    #   ```
    # Which simply results in:
    start_idcs -= 1

    ranges = []
    for interpro_accession, start_idx, end_idx in zip(
        interpro_accessions, start_idcs, end_idcs
    ):
        # NOTE: Skip unintegrated Interpro labels, for now.
        if interpro_accession == "-":
            continue

        ranges.append(
            InterProRangeAnnotation(
                interpro_accession=interpro_accession,
                start_idx=start_idx,
                end_idx=end_idx,
            )
        )

    return ranges
