from __future__ import annotations

from collections.abc import Sequence

import pandas as pd
from localcider.sequenceParameters import SequenceParameters

BASE_FEATURE_COLUMNS = [
    "Hydropathy",
    "FCR",
    "IDR",
    "LCR",
    "PScore",
    "PLAAC",
    "catGRANULE",
    "DeepCoil",
]

HUMAN_FEATURE_COLUMNS = BASE_FEATURE_COLUMNS + ["Phos freq", "DeepPhase"]

# DeepCoil v2022 definition (FEATURES.md §8): binary indicator with
# threshold 0.82 on the raw coiled-coil probability.
DEEPCOIL_THRESHOLD = 0.82

# catGRANULE residue profile window (CATGRANULE_SPEC §5 / FEATURES.md §7):
# a 51-residue centered window, covering the interior positions
# 25..L-26 (1-indexed) → length L-50, offset 24.
CATGRANULE_RESIDUE_WINDOW = 51
CATGRANULE_RESIDUE_OFFSET = 24

KYTE_DOOLITTLE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

CHARGED_RESIDUES = frozenset({"D", "E", "K", "R"})


class FeatureSchemaError(ValueError):
    """Raised when feature inputs do not satisfy a model schema."""


def compute_native_features(sequence: str) -> dict[str, float]:
    normalized = sequence.upper()
    if not normalized:
        raise FeatureSchemaError("Cannot compute features for an empty sequence")
    unknown = sorted({residue for residue in normalized if residue not in KYTE_DOOLITTLE})
    if unknown:
        raise FeatureSchemaError(f"Unsupported residues in sequence: {', '.join(unknown)}")
    length = len(normalized)
    return {
        "length": float(length),
        "Hydropathy": SequenceParameters(normalized).get_uversky_hydropathy(),
        "FCR": sum(1 for residue in normalized if residue in CHARGED_RESIDUES) / length,
    }


def build_feature_matrix(
    sequence_records: pd.DataFrame,
    imported_features: pd.DataFrame | None = None,
    *,
    feature_columns: Sequence[str] = BASE_FEATURE_COLUMNS,
) -> pd.DataFrame:
    _require_columns(sequence_records, {"UniprotEntry", "sequence"}, "sequence records")
    native_rows = []
    for row in sequence_records.itertuples(index=False):
        features = compute_native_features(str(row.sequence))
        native_rows.append({"UniprotEntry": str(row.UniprotEntry), **features})
    native = pd.DataFrame(native_rows)

    if imported_features is None:
        matrix = native
    else:
        _require_columns(imported_features, {"UniprotEntry"}, "imported features")
        matrix = native.merge(
            imported_features,
            on="UniprotEntry",
            how="left",
            suffixes=("", "_imported"),
        )
        for column in ["Hydropathy", "FCR", "length"]:
            imported_column = f"{column}_imported"
            if imported_column in matrix.columns:
                matrix = matrix.drop(columns=[imported_column])

    missing = [column for column in feature_columns if column not in matrix.columns]
    if missing:
        raise FeatureSchemaError(f"Missing required feature columns: {', '.join(missing)}")
    return matrix


def validate_feature_columns(frame: pd.DataFrame, feature_columns: Sequence[str]) -> None:
    missing = [column for column in feature_columns if column not in frame.columns]
    if missing:
        raise FeatureSchemaError(f"Missing required feature columns: {', '.join(missing)}")


def _require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise FeatureSchemaError(f"{name} missing columns: {', '.join(missing)}")


# ---------------------------------------------------------------------------
# Residue-level alignment utilities (v2022 definitions)
# ---------------------------------------------------------------------------


def catgranule_profile_alignment(sequence_length: int) -> range:
    """Position range (0-indexed) covered by a catGRANULE residue profile.

    The profile is a 51-residue sliding-window mean (radius 25) over the
    interior span: 0-indexed ``24 .. L-27`` (i.e. offset 24, length ``L-50``;
    FEATURES.md §7 / CATGRANULE_SPEC §5). Empty range when the sequence is
    too short for the window.
    """
    radius = (CATGRANULE_RESIDUE_WINDOW - 1) // 2
    start = radius - 1
    stop = sequence_length - radius - 1
    if stop <= start:
        return range(0, 0)
    return range(start, stop)


def pscore_alignment(sequence_length: int, residue_length: int) -> range:
    """Position range (0-indexed) covered by a PScore residue array.

    PScore omits a sequence-dependent number of head residues (empirically
    1..L-delta; FEATURES.md §5): align the shorter array to the *suffix* of
    the full sequence. Returns an empty range when the residue array is
    longer than the sequence.
    """
    if residue_length < 0 or residue_length > sequence_length:
        return range(0, 0)
    return range(sequence_length - residue_length, sequence_length)


def binary_deepcoil(raw_cc: Sequence[float]) -> float:
    """v2022 DeepCoil feature value: ``1.0 if max(raw_cc) >= 0.82 else 0.0``."""
    return 1.0 if max(raw_cc) >= DEEPCOIL_THRESHOLD else 0.0
