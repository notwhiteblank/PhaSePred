"""Bit-exact fidelity of data/processed/*.tsv to the paper workbooks (E1 G1-G3).

Needs the private supplementary workbooks, so the whole module skips with an
explicit reason when they are absent (public checkout, CI). The checks that do
not need the workbooks -- sizes, digests, schema, duplicates, manifest -- live in
tests/test_training_data_integrity.py and never skip.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phasepred.paper import load_paper_features
from phasepred.training_data import (
    COLUMNS,
    EXPECTED_ROWS,
    FLOAT_COLUMNS,
    INT_COLUMNS,
    SCOPES,
    SPLITS,
    TASKS_BY_SCOPE,
    read_manifest,
    read_table,
    sha256_file,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
S2 = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd02.xlsx"
S3 = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd03.xlsx"
KEYS = [(scope, split) for scope in SCOPES for split in SPLITS]

pytestmark = pytest.mark.skipif(
    not (S2.is_file() and S3.is_file()),
    reason="paper supplementary workbooks absent (see docs/DATA_SOURCES.md)",
)


@pytest.fixture(scope="module")
def paper_features() -> dict[str, pd.DataFrame]:
    features = load_paper_features(S2, S3)
    return {"train": features.train, "test": features.test}


def scope_slice(paper_features: dict[str, pd.DataFrame], scope: str, split: str):
    source = paper_features[split]
    return source[source["task"].astype(str).isin(TASKS_BY_SCOPE[scope])]


def test_row_counts_and_order_match_the_paper_sheets(paper_features):
    for scope, split in KEYS:
        want = scope_slice(paper_features, scope, split)
        got = read_table(scope, split)
        assert len(got) == len(want) == EXPECTED_ROWS[(scope, split)]
        assert got["UniprotEntry"].tolist() == want["UniprotEntry"].astype(str).tolist()


def test_values_are_bit_identical(paper_features):
    for scope, split in KEYS:
        want = scope_slice(paper_features, scope, split).reset_index(drop=True)
        got = read_table(scope, split)
        for column in COLUMNS[(scope, split)]:
            where = (scope, split, column)
            if column in ("task", "split", "label", "source_sheet"):
                got_text = got[column].astype(str).tolist()
                want_text = want[column].astype(str).tolist()
                assert got_text == want_text, where
            elif column in INT_COLUMNS or column in FLOAT_COLUMNS:
                gv = pd.to_numeric(got[column], errors="coerce").to_numpy(dtype="float64")
                wv = pd.to_numeric(want[column], errors="coerce").to_numpy(dtype="float64")
                assert np.array_equal(gv, wv, equal_nan=True), where
                if column in FLOAT_COLUMNS:
                    assert gv.tobytes() == wv.tobytes(), where
            else:
                w = want[column].astype("string")
                w = w.mask(w == "", pd.NA)
                assert got[column].isna().tolist() == w.isna().tolist(), where
                assert got[column].dropna().tolist() == w.dropna().tolist(), where


def test_no_imputation_happened(paper_features):
    for scope, split in KEYS:
        want = scope_slice(paper_features, scope, split)
        got = read_table(scope, split)
        for column in COLUMNS[(scope, split)]:
            if column in ("task", "split", "label", "source_sheet") or column not in want:
                continue
            source = want[column]
            if column in INT_COLUMNS or column in FLOAT_COLUMNS:
                na_want = int(pd.to_numeric(source, errors="coerce").isna().sum())
            else:
                text = source.astype("string")
                na_want = int((text.isna() | (text == "")).sum())
            assert int(got[column].isna().sum()) == na_want, (scope, split, column)


def test_manifest_records_the_real_workbook_digests():
    manifest = read_manifest().set_index("file")
    want_s2 = sha256_file(S2)
    want_s3 = sha256_file(S3)
    for name, row in manifest.iterrows():
        assert row["sd02_sha256"] == want_s2, name
        assert row["sd03_sha256"] == want_s3, name
