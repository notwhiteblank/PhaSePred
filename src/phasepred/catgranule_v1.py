from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

import pandas as pd

CATGRANULE_WEIGHTS = MappingProxyType(
    {
        "rc": 0.48,
        "rn": 7.24,
        "dc": 0.26,
        "dn": 11.54,
        "prg": 1.98,
        "pfg": 1.42,
        "length": 0.25,
    }
)
WINDOW_RADIUS = 3


class CatGranuleInputError(ValueError):
    """Raised when a sequence cannot be scored by the reconstruction."""


@dataclass(frozen=True)
class CatGranuleScale:
    rc: float
    rn: float
    dc: float
    dn: float
    pfg: float
    prg: float


@dataclass(frozen=True)
class CatGranuleNormalization:
    mean: float
    standard_deviation: float

    def __post_init__(self) -> None:
        if self.standard_deviation <= 0:
            raise CatGranuleInputError("Normalization standard deviation must be positive")


TABLE_S4_SCALES: Mapping[str, CatGranuleScale] = MappingProxyType(
    {
        "A": CatGranuleScale(rc=0.45, rn=0.46, dc=0.08, dn=0.03, pfg=0.00, prg=0.00),
        "C": CatGranuleScale(rc=0.00, rn=0.00, dc=0.24, dn=0.01, pfg=0.00, prg=0.00),
        "D": CatGranuleScale(rc=0.61, rn=0.61, dc=0.56, dn=0.29, pfg=0.00, prg=0.00),
        "E": CatGranuleScale(rc=0.44, rn=0.70, dc=0.00, dn=0.09, pfg=0.00, prg=0.00),
        "F": CatGranuleScale(rc=0.59, rn=0.38, dc=0.05, dn=0.08, pfg=1.00, prg=0.00),
        "G": CatGranuleScale(rc=0.92, rn=0.52, dc=0.63, dn=1.00, pfg=1.00, prg=1.00),
        "H": CatGranuleScale(rc=0.51, rn=0.38, dc=0.40, dn=0.19, pfg=0.00, prg=0.00),
        "I": CatGranuleScale(rc=0.46, rn=0.37, dc=0.17, dn=0.09, pfg=0.00, prg=0.00),
        "K": CatGranuleScale(rc=0.59, rn=1.00, dc=0.18, dn=0.22, pfg=0.00, prg=0.00),
        "L": CatGranuleScale(rc=0.30, rn=0.37, dc=0.06, dn=0.05, pfg=0.00, prg=0.00),
        "M": CatGranuleScale(rc=0.24, rn=0.05, dc=0.06, dn=0.11, pfg=0.00, prg=0.00),
        "N": CatGranuleScale(rc=0.67, rn=0.40, dc=0.53, dn=0.67, pfg=0.00, prg=0.00),
        "P": CatGranuleScale(rc=0.69, rn=0.40, dc=1.00, dn=0.09, pfg=0.00, prg=0.00),
        "Q": CatGranuleScale(rc=0.56, rn=0.51, dc=0.24, dn=0.21, pfg=0.00, prg=0.00),
        "R": CatGranuleScale(rc=0.89, rn=0.72, dc=0.17, dn=0.05, pfg=0.00, prg=1.00),
        "S": CatGranuleScale(rc=0.70, rn=0.44, dc=0.47, dn=0.24, pfg=0.00, prg=0.00),
        "T": CatGranuleScale(rc=0.44, rn=0.45, dc=0.50, dn=0.00, pfg=0.00, prg=0.00),
        "V": CatGranuleScale(rc=0.44, rn=0.39, dc=0.01, dn=0.14, pfg=0.00, prg=0.00),
        "W": CatGranuleScale(rc=0.46, rn=0.21, dc=0.23, dn=0.09, pfg=0.00, prg=0.00),
        "Y": CatGranuleScale(rc=1.00, rn=0.54, dc=0.45, dn=0.17, pfg=0.00, prg=0.00),
    }
)


TABLE_S4_SCALE_COLUMNS = ("rc", "rn", "dc", "dn", "pfg", "prg")

TABLE_S1_VALIDATION_SETS = frozenset({"Granule Forming", "Granule Related"})
_TABLE_S1_YEAST_ID_RE = re.compile(r"^(?P<uniprot_entry>[A-Z0-9]+_YEAST)(?:\s+(?P<source>.+))?$")


def extract_table_s1_labels(input_path: Path) -> pd.DataFrame:
    raw = pd.read_excel(input_path, engine="calamine", sheet_name=0, header=None)
    rows: list[dict[str, str]] = []

    for row_number, row in raw.iloc[1:].iterrows():
        gene_cell = "" if pd.isna(row.iloc[0]) else str(row.iloc[0]).strip()
        set_name = "" if pd.isna(row.iloc[2]) else str(row.iloc[2]).strip()
        if not gene_cell and not set_name:
            continue
        if set_name not in TABLE_S1_VALIDATION_SETS:
            raise ValueError(
                f"Unexpected Table S1 set label at Excel row {row_number + 1}: {set_name!r}"
            )

        match = _TABLE_S1_YEAST_ID_RE.match(gene_cell)
        if match is None:
            raise ValueError(
                f"Could not parse Table S1 gene cell at Excel row {row_number + 1}: {gene_cell!r}"
            )

        rows.append(
            {
                "uniprot_entry_name": match.group("uniprot_entry"),
                "table_s1_set": set_name,
                "source_note": match.group("source") or "",
            }
        )

    labels = pd.DataFrame(rows, columns=["uniprot_entry_name", "table_s1_set", "source_note"])
    duplicates = labels["uniprot_entry_name"][labels["uniprot_entry_name"].duplicated()].unique()
    if len(duplicates):
        raise ValueError(f"Duplicate Table S1 entry names: {', '.join(sorted(duplicates))}")
    return labels


def centered_window(sequence: str, index: int, *, radius: int = WINDOW_RADIUS) -> str:
    normalized = _normalize_sequence(sequence)
    if index < 0 or index >= len(normalized):
        raise CatGranuleInputError(f"Sequence index out of range: {index}")
    start = max(0, index - radius)
    stop = min(len(normalized), index + radius + 1)
    return normalized[start:stop]


def residue_score(sequence: str, index: int) -> float:
    window = centered_window(sequence, index)
    window_scales = [TABLE_S4_SCALES[residue] for residue in window]
    width = len(window_scales)
    means = {
        column: sum(getattr(scale, column) for scale in window_scales) / width
        for column in TABLE_S4_SCALE_COLUMNS
    }
    return (
        CATGRANULE_WEIGHTS["rc"] * means["rc"]
        + CATGRANULE_WEIGHTS["rn"] * means["rn"]
        + CATGRANULE_WEIGHTS["dc"] * means["dc"]
        + CATGRANULE_WEIGHTS["dn"] * means["dn"]
        + CATGRANULE_WEIGHTS["prg"] * means["prg"]
        + CATGRANULE_WEIGHTS["pfg"] * means["pfg"]
    )


def residue_scores(sequence: str) -> list[float]:
    normalized = _normalize_sequence(sequence)
    return [residue_score(normalized, index) for index in range(len(normalized))]


def raw_score(sequence: str, *, include_length_term: bool = True) -> float:
    normalized = _normalize_sequence(sequence)
    scores = residue_scores(normalized)
    score = sum(scores) / len(scores)
    if include_length_term:
        score += CATGRANULE_WEIGHTS["length"] * math.log(len(normalized))
    return score


def normalize_score(score: float, normalization: CatGranuleNormalization) -> float:
    return (score - normalization.mean) / normalization.standard_deviation


def profile(
    sequence: str,
    *,
    window_size: int = 50,
    normalization: CatGranuleNormalization | None = None,
) -> list[float]:
    if window_size <= 0:
        raise CatGranuleInputError("Profile window size must be positive")

    scores = residue_scores(sequence)
    radius_left = (window_size - 1) // 2
    radius_right = window_size // 2
    values: list[float] = []
    for index in range(len(scores)):
        start = max(0, index - radius_left)
        stop = min(len(scores), index + radius_right + 1)
        value = sum(scores[start:stop]) / (stop - start)
        if normalization is not None:
            value = normalize_score(value, normalization)
        values.append(value)
    return values


def granule_strength(normalized_profile: list[float]) -> float:
    positive_values = [value for value in normalized_profile if value > 0]
    if not positive_values:
        return 0.0
    return sum(positive_values) / len(positive_values)


def build_entry_name_index(
    cache_path: str | Path,
) -> dict[str, str]:
    """Build a reverse index from UniProt entry name (e.g. BEM2_YEAST) to accession.

    Reads a JSONL sequence cache and parses the UniProt header of each record.
    Only maps entries where the header conforms to the standard ``>sp|ACC|NAME``
    or ``>tr|ACC|NAME`` format.
    """
    from phasepred.sequences import SequenceCache

    cache = SequenceCache(cache_path)
    index: dict[str, str] = {}
    for accession, record in cache._read_all().items():
        parts = record.header[1:].split(maxsplit=1)[0].split("|")
        if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
            index[parts[2]] = accession
    return index


def _normalize_sequence(sequence: str) -> str:
    normalized = sequence.upper()
    if not normalized:
        raise CatGranuleInputError("Cannot score an empty sequence")
    unknown = sorted({residue for residue in normalized if residue not in TABLE_S4_SCALES})
    if unknown:
        raise CatGranuleInputError(f"Unsupported residues in sequence: {', '.join(unknown)}")
    return normalized
