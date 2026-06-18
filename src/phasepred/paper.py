from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


class PaperSplitError(ValueError):
    """Raised when paper supplemental split tables cannot be parsed."""


@dataclass(frozen=True)
class PaperSplit:
    train: pd.DataFrame
    test: pd.DataFrame

    @property
    def all_rows(self) -> pd.DataFrame:
        return pd.concat([self.train, self.test], ignore_index=True)


@dataclass(frozen=True)
class PaperFeatures:
    train: pd.DataFrame
    test: pd.DataFrame

    @property
    def all_rows(self) -> pd.DataFrame:
        return pd.concat([self.train, self.test], ignore_index=True)


TRAIN_SHEETS: dict[str, tuple[str, str]] = {
    "SaPS": ("SaPS", "NoPS"),
    "PdPS": ("PdPS", "NoPS"),
    "hSaPS": ("hSaPS", "hNoPS"),
    "hPdPS": ("hPdPS", "hNoPS"),
}

TEST_SHEETS: dict[str, tuple[tuple[str, ...], str]] = {
    "SaPS": (("SaPS-test", "PS-test"), "NoPS-test"),
    "PdPS": (("PdPS-test", "PS-test"), "NoPS-test"),
    "hSaPS": (("hSaPS-test", "hPS-test"), "hNoPS-test"),
    "hPdPS": (("hPdPS-test", "hPS-test"), "hNoPS-test"),
}


def load_paper_split(s2_path: str | Path, s3_path: str | Path) -> PaperSplit:
    s2 = Path(s2_path)
    s3 = Path(s3_path)
    if not s2.exists():
        raise PaperSplitError(f"Dataset S2 workbook not found: {s2}")
    if not s3.exists():
        raise PaperSplitError(f"Dataset S3 workbook not found: {s3}")

    train_rows: list[dict[str, object]] = []
    for task, (positive_sheet, negative_sheet) in TRAIN_SHEETS.items():
        train_rows.extend(_rows_for_sheet(s2, positive_sheet, task, "train", 1))
        train_rows.extend(_rows_for_sheet(s2, negative_sheet, task, "train", 0))

    test_rows: list[dict[str, object]] = []
    for task, (positive_sheets, negative_sheet) in TEST_SHEETS.items():
        for positive_sheet in positive_sheets:
            test_rows.extend(_rows_for_sheet(s3, positive_sheet, task, "test", 1))
        test_rows.extend(_rows_for_sheet(s3, negative_sheet, task, "test", 0))

    return PaperSplit(
        train=pd.DataFrame(train_rows, columns=_split_columns()),
        test=pd.DataFrame(test_rows, columns=_split_columns()),
    )


def load_paper_features(s2_path: str | Path, s3_path: str | Path) -> PaperFeatures:
    s2 = Path(s2_path)
    s3 = Path(s3_path)
    if not s2.exists():
        raise PaperSplitError(f"Dataset S2 workbook not found: {s2}")
    if not s3.exists():
        raise PaperSplitError(f"Dataset S3 workbook not found: {s3}")

    train_rows: list[dict[str, object]] = []
    for task, (positive_sheet, negative_sheet) in TRAIN_SHEETS.items():
        train_rows.extend(_feature_rows_for_sheet(s2, positive_sheet, task, "train", 1))
        train_rows.extend(_feature_rows_for_sheet(s2, negative_sheet, task, "train", 0))

    test_rows: list[dict[str, object]] = []
    for task, (positive_sheets, negative_sheet) in TEST_SHEETS.items():
        for positive_sheet in positive_sheets:
            test_rows.extend(_feature_rows_for_sheet(s3, positive_sheet, task, "test", 1))
        test_rows.extend(_feature_rows_for_sheet(s3, negative_sheet, task, "test", 0))

    return PaperFeatures(
        train=pd.DataFrame(train_rows),
        test=pd.DataFrame(test_rows),
    )


def _split_columns() -> list[str]:
    return ["task", "split", "label", "UniprotEntry", "source_sheet"]


def _rows_for_sheet(
    workbook: Path,
    sheet_name: str,
    task: str,
    split: str,
    label: int,
) -> list[dict[str, object]]:
    frame = _read_sheet(workbook, sheet_name)
    return [
        {
            "task": task,
            "split": split,
            "label": label,
            "UniprotEntry": accession,
            "source_sheet": sheet_name,
        }
        for accession in (frame["UniprotEntry"].dropna().astype(str).str.strip().tolist())
        if accession
    ]


def _feature_rows_for_sheet(
    workbook: Path,
    sheet_name: str,
    task: str,
    split: str,
    label: int,
) -> list[dict[str, object]]:
    frame = _read_sheet(workbook, sheet_name)
    records = frame.to_dict(orient="records")
    rows: list[dict[str, object]] = []
    for record in records:
        accession = str(record.get("UniprotEntry") or "").strip()
        if not accession:
            continue
        row = dict(record)
        row.update(
            {
                "task": task,
                "split": split,
                "label": label,
                "UniprotEntry": accession,
                "source_sheet": sheet_name,
            }
        )
        rows.append(row)
    return rows


def _read_sheet(workbook: Path, sheet_name: str) -> pd.DataFrame:
    try:
        frame = pd.read_excel(workbook, sheet_name=sheet_name)
    except ValueError as exc:
        raise PaperSplitError(f"Required sheet {sheet_name!r} not found in {workbook}") from exc
    if "UniprotEntry" not in frame.columns:
        frame = pd.read_excel(workbook, sheet_name=sheet_name, header=1)
    if "UniprotEntry" not in frame.columns:
        raise PaperSplitError(f"Sheet {sheet_name!r} in {workbook} lacks UniprotEntry column")
    return frame
