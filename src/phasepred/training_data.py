"""On-disk contract for the PhaSePred training tables.

This module owns everything about ``data/processed/*.tsv``: the column schema,
the dtypes, the filenames, the audited sizes, and both directions of the
serialization. ``scripts/build_training_tsv.py`` writes the tables through
:func:`select_and_cast` and :func:`write_tsv`; ``scripts/train_phasepred.py``
reads them through :func:`read_table` and :func:`load_task_frame`.

The tables are a faithful serialization of
:func:`phasepred.paper.load_paper_features`. Row order and row multiplicity are
preserved exactly, missing values stay missing (the paper protocol passes NaN
straight to XGBoost), and no imputation, sorting or de-duplication is applied.
Negatives repeat across tasks because the paper's ``NoPS`` / ``hNoPS`` sheets
serve two modes each; that multiplicity is part of the protocol.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from phasepred.data import data_path
from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS

SCOPES = ("base", "human")
SPLITS = ("train", "test")
TASKS_BY_SCOPE: dict[str, tuple[str, ...]] = {
    "base": ("SaPS", "PdPS"),
    "human": ("hSaPS", "hPdPS"),
}

#: Bump when a serialization change alters the output bytes.
SCRIPT_VERSION = "1.0.0"
#: When the supplementary workbooks were obtained. A constant rather than a
#: clock reading: nothing wall-clock may enter the output, or repeated builds
#: would differ and ``--check`` would drift.
DATA_ACCESSION_DATE = "2026-05-06"
SOURCE_WORKBOOKS = (
    "pnas.2115369119.sd02.xlsx",
    "pnas.2115369119.sd03.xlsx",
)

META_COLUMNS = ["UniprotEntry", "task", "split", "label", "source_sheet"]
PAPER_COLUMNS = ["Gene name", "Organism", "Organism ID", "length"]
TEST_ONLY_COLUMNS = ["Source"]

INT_COLUMNS = frozenset({"label", "Organism ID", "length", "DeepCoil"})
FLOAT_COLUMNS = (
    frozenset(BASE_FEATURE_COLUMNS) | frozenset(HUMAN_FEATURE_COLUMNS)
) - INT_COLUMNS
TEXT_COLUMNS = frozenset(
    {"UniprotEntry", "task", "split", "source_sheet", "Gene name", "Organism", "Source"}
)

_BASE_TRAIN = META_COLUMNS + PAPER_COLUMNS + list(BASE_FEATURE_COLUMNS)
_HUMAN_TRAIN = META_COLUMNS + PAPER_COLUMNS + list(HUMAN_FEATURE_COLUMNS)
COLUMNS: dict[tuple[str, str], list[str]] = {
    ("base", "train"): _BASE_TRAIN,
    ("human", "train"): _HUMAN_TRAIN,
    ("base", "test"): TEST_ONLY_COLUMNS + _BASE_TRAIN,
    ("human", "test"): TEST_ONLY_COLUMNS + _HUMAN_TRAIN,
}

FILENAMES: dict[tuple[str, str], str] = {
    ("base", "train"): "chen2022_s2s3_base-features_train.tsv",
    ("base", "test"): "chen2022_s2s3_base-features_test.tsv",
    ("human", "train"): "chen2022_s2s3_human-features_train.tsv",
    ("human", "test"): "chen2022_s2s3_human-features_test.tsv",
}

#: Audited 2026-09-16 against the supplementary workbooks; PLAN-E2E §5 E1.
EXPECTED_ROWS: dict[tuple[str, str], int] = {
    ("base", "train"): 96658,
    ("base", "test"): 24416,
    ("human", "train"): 17757,
    ("human", "test"): 4540,
}
EXPECTED_BYTES: dict[tuple[str, str], int] = {
    ("base", "train"): 13399739,
    ("base", "test"): 3508476,
    ("human", "train"): 2800875,
    ("human", "test"): 736518,
}
#: Frozen from the audited generation of 2026-09-16; the workbook digests are
#: recorded in data/processed/MANIFEST.tsv. Any change to the tables must be
#: deliberate: bump SCRIPT_VERSION, regenerate, re-freeze these digests, and
#: update EXPECTED_BYTES.
EXPECTED_SHA256: dict[tuple[str, str], str] = {
    ("base", "train"): "df9081b001936adb678de36ac186811df7ed0f5db1e09079e6e788049d385bb7",
    ("base", "test"): "b26253f921217375f4873d5714058fcf20925729b21b8381a4eba27f3940ff9a",
    ("human", "train"): "fedecb3850263de9b561a6a5685f91b3f1235135b4baf6aca4504031969cad59",
    ("human", "test"): "cb36429d16f575785213dc625708c302e093129e4dc29ddad60facad36928134",
}

MANIFEST_FILENAME = "MANIFEST.tsv"
MANIFEST_COLUMNS = [
    "file",
    "scope",
    "split",
    "rows",
    "columns",
    "bytes",
    "sha256",
    "sd02_sha256",
    "sd03_sha256",
    "script",
    "script_version",
    "data_accession_date",
]
BUILD_SCRIPT = "scripts/build_training_tsv.py"


class TrainingDataSchemaError(ValueError):
    """Raised when a training table violates the documented schema."""


def _key(scope: str, split: str) -> tuple[str, str]:
    if scope not in SCOPES:
        raise TrainingDataSchemaError(f"Unknown scope {scope!r}; expected one of {SCOPES}")
    if split not in SPLITS:
        raise TrainingDataSchemaError(f"Unknown split {split!r}; expected one of {SPLITS}")
    return (scope, split)


def processed_dir(root: Path | None = None) -> Path:
    """Directory holding the four tables.

    ``root`` overrides the directory itself (tests, ``--output-dir``). Without
    it the tables live at ``<data root>/data/processed``.
    """
    if root is not None:
        return Path(root)
    return data_path("data", "processed")


def table_path(scope: str, split: str, *, root: Path | None = None) -> Path:
    return processed_dir(root) / FILENAMES[_key(scope, split)]


def manifest_path(*, root: Path | None = None) -> Path:
    return processed_dir(root) / MANIFEST_FILENAME


def dtype_map(scope: str, split: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for column in COLUMNS[_key(scope, split)]:
        if column in INT_COLUMNS:
            mapping[column] = "Int64"
        elif column in FLOAT_COLUMNS:
            mapping[column] = "float64"
        else:
            mapping[column] = "string"
    return mapping


def validate_schema(frame: pd.DataFrame, scope: str, split: str) -> None:
    """Raise :class:`TrainingDataSchemaError` unless ``frame`` fits the contract."""
    key = _key(scope, split)
    name = FILENAMES[key]
    expected = COLUMNS[key]
    if list(frame.columns) != expected:
        raise TrainingDataSchemaError(
            f"{name}: column mismatch.\n"
            f"  expected: {expected}\n"
            f"  got:      {list(frame.columns)}"
        )
    allowed = set(TASKS_BY_SCOPE[scope])
    got_tasks = set(frame["task"].astype(str))
    if not got_tasks <= allowed:
        raise TrainingDataSchemaError(f"{name}: out-of-scope tasks {sorted(got_tasks - allowed)}")
    got_splits = set(frame["split"].astype(str))
    if got_splits != {split}:
        raise TrainingDataSchemaError(
            f"{name}: expected split={split!r} only, got {sorted(got_splits)}"
        )
    labels = set(pd.to_numeric(frame["label"], errors="coerce").dropna().unique().tolist())
    if not labels <= {0, 1}:
        raise TrainingDataSchemaError(f"{name}: label must be 0/1, got {sorted(labels)}")
    duplicated = int(frame.duplicated(subset=["task", "split", "UniprotEntry"]).sum())
    if duplicated:
        raise TrainingDataSchemaError(
            f"{name}: {duplicated} duplicated (task, split, UniprotEntry) rows. "
            "The paper sheets contain none, so the source workbooks changed."
        )


def select_and_cast(frame: pd.DataFrame, scope: str, split: str) -> pd.DataFrame:
    """Project ``load_paper_features()`` output onto one table's schema.

    Selects by exact column name (the ``NoPS`` / ``hNoPS`` sheets carry unnamed
    extra columns), casts to the declared dtypes, normalizes empty strings in
    text columns to missing, and preserves row order and multiplicity.
    """
    key = _key(scope, split)
    columns = COLUMNS[key]
    keep = frame[frame["task"].astype(str).isin(TASKS_BY_SCOPE[scope])]
    keep = keep[keep["split"].astype(str) == split]
    if keep.empty:
        return keep.reindex(columns=columns)
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise TrainingDataSchemaError(
            f"source frame lacks required columns: {missing} (scope={scope}, split={split})"
        )
    out = keep.loc[:, columns].copy()
    for column in columns:
        if column in INT_COLUMNS:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
        elif column in FLOAT_COLUMNS:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("float64")
        else:
            text = out[column].astype("string")
            out[column] = text.mask(text == "", pd.NA)
    return out


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Write ``frame`` as TSV: UTF-8, LF, header, missing values as empty fields."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path, sep="\t", index=False, na_rep="", lineterminator="\n", encoding="utf-8"
    )


def read_tsv(path: Path, scope: str, split: str) -> pd.DataFrame:
    """Read a table back and validate it.

    ``keep_default_na=False`` with an explicit ``na_values=[""]`` matters: only
    a truly empty field counts as missing, so a gene symbol that literally reads
    "NA" survives the round trip.

    ``float_precision="round_trip"`` matters just as much: the workbooks carry
    17-significant-digit floats and pandas' default parser ("high") is off by
    one ULP on them, which would break the bit-exactness contract.
    """
    frame = pd.read_csv(
        Path(path),
        sep="\t",
        dtype=dtype_map(scope, split),
        keep_default_na=False,
        na_values=[""],
        float_precision="round_trip",
        encoding="utf-8",
    )
    validate_schema(frame, scope, split)
    return frame


def read_table(scope: str, split: str, *, root: Path | None = None) -> pd.DataFrame:
    path = table_path(scope, split, root=root)
    if not path.is_file():
        raise TrainingDataSchemaError(
            f"Training table missing: {path}. Build it with "
            f"`python {BUILD_SCRIPT}` inside a checkout that has the paper "
            "supplementary workbooks (see docs/DATA_SOURCES.md)."
        )
    return read_tsv(path, scope, split)


def load_task_frame(
    scope: str, split: str, task: str, *, root: Path | None = None
) -> pd.DataFrame:
    """One task's rows, in source order, with a reset positional index."""
    _key(scope, split)
    if task not in TASKS_BY_SCOPE[scope]:
        raise TrainingDataSchemaError(
            f"task {task!r} does not belong to scope {scope!r} "
            f"(expected one of {TASKS_BY_SCOPE[scope]})"
        )
    frame = read_table(scope, split, root=root)
    return frame[frame["task"].astype(str) == task].reset_index(drop=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_expected_size(scope: str, split: str, *, root: Path | None = None) -> None:
    """Assert bytes, digest and row count against the audited values.

    Checked in that order so a truncated or appended file is reported as a size
    drift before the more expensive parse.
    """
    key = _key(scope, split)
    path = table_path(scope, split, root=root)
    if not path.is_file():
        raise TrainingDataSchemaError(f"Training table missing: {path}")
    size = path.stat().st_size
    if size != EXPECTED_BYTES[key]:
        raise TrainingDataSchemaError(
            f"{path.name}: {size} bytes, expected {EXPECTED_BYTES[key]}"
        )
    expected_sha = EXPECTED_SHA256.get(key)
    if expected_sha:
        actual_sha = sha256_file(path)
        if actual_sha != expected_sha:
            raise TrainingDataSchemaError(
                f"{path.name}: sha256 {actual_sha}, expected {expected_sha}"
            )
    rows = len(read_tsv(path, scope, split))
    if rows != EXPECTED_ROWS[key]:
        raise TrainingDataSchemaError(
            f"{path.name}: {rows} rows, expected {EXPECTED_ROWS[key]}"
        )


def manifest_row(
    path: Path,
    scope: str,
    split: str,
    *,
    rows: int,
    sd02_sha256: str,
    sd03_sha256: str,
) -> dict[str, object]:
    path = Path(path)
    key = _key(scope, split)
    return {
        "file": path.name,
        "scope": scope,
        "split": split,
        "rows": rows,
        "columns": len(COLUMNS[key]),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "sd02_sha256": sd02_sha256,
        "sd03_sha256": sd03_sha256,
        "script": BUILD_SCRIPT,
        "script_version": SCRIPT_VERSION,
        "data_accession_date": DATA_ACCESSION_DATE,
    }


def write_manifest(rows: Sequence[dict[str, object]], path: Path) -> None:
    frame = pd.DataFrame(list(rows), columns=MANIFEST_COLUMNS)
    write_tsv(frame, Path(path))


def read_manifest(*, root: Path | None = None) -> pd.DataFrame:
    path = manifest_path(root=root)
    if not path.is_file():
        raise TrainingDataSchemaError(f"Manifest missing: {path}")
    return pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False, na_values=[""], encoding="utf-8"
    )
