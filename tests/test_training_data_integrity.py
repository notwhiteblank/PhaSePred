"""Integrity checks for the committed data/processed tables (E1).

The four tables, the manifest and the frozen digests are committed artifacts, so
these checks must run everywhere -- a public checkout, CI -- with no private
workbooks and with no skip of any kind. A missing or drifted table fails loudly
here rather than skipping quietly.

The workbook-dependent fidelity checks (row order, bit-exact values, no
imputation) and the real workbook-digest check live in
tests/test_training_data_real.py, which is gated on the private supplementary
workbooks.
"""

from __future__ import annotations

from phasepred.training_data import (
    COLUMNS,
    EXPECTED_BYTES,
    EXPECTED_ROWS,
    EXPECTED_SHA256,
    FILENAMES,
    SCOPES,
    SPLITS,
    read_manifest,
    read_table,
    sha256_file,
    table_path,
    validate_expected_size,
)

KEYS = [(scope, split) for scope in SCOPES for split in SPLITS]

#: Audited 2026-09-16: duplicates come only from sheets shared between two
#: tasks (NoPS / hNoPS / PS-test / NoPS-test / hNoPS-test).
EXPECTED_DUPLICATES = {
    ("base", "train"): 48158,
    ("base", "test"): 12115,
    ("human", "train"): 8801,
    ("human", "test"): 2223,
}


def test_tables_exist_and_match_the_audited_sizes():
    for scope, split in KEYS:
        assert table_path(scope, split).is_file(), FILENAMES[(scope, split)]
        validate_expected_size(scope, split)


def test_sha256_is_frozen_and_matches():
    assert set(EXPECTED_SHA256) == set(KEYS)
    for scope, split in KEYS:
        assert sha256_file(table_path(scope, split)) == EXPECTED_SHA256[(scope, split)]


def test_missing_values_are_empty_fields_not_sentinels():
    for scope, split in KEYS:
        text = table_path(scope, split).read_text(encoding="utf-8")
        for sentinel in ("NaN", "<NA>", "NULL", "None", "nan"):
            assert f"\t{sentinel}\t" not in text, (FILENAMES[(scope, split)], sentinel)
            assert f"\n{sentinel}\t" not in text, (FILENAMES[(scope, split)], sentinel)
            assert f"\t{sentinel}\n" not in text, (FILENAMES[(scope, split)], sentinel)


def test_no_duplicate_accession_within_a_task_and_split():
    for scope, split in KEYS:
        got = read_table(scope, split)
        assert not got.duplicated(subset=["task", "split", "UniprotEntry"]).any()


def test_duplicate_rows_are_cross_task_only():
    for key, count in EXPECTED_DUPLICATES.items():
        got = read_table(*key)
        actual = len(got) - got["UniprotEntry"].nunique()
        assert actual == count, f"{key}: {actual} duplicates, expected {count}"


def test_manifest_matches_the_files_on_disk():
    manifest = read_manifest().set_index("file")
    assert len(manifest) == 4
    for scope, split in KEYS:
        name = FILENAMES[(scope, split)]
        path = table_path(scope, split)
        row = manifest.loc[name]
        assert int(row["rows"]) == EXPECTED_ROWS[(scope, split)]
        assert int(row["bytes"]) == path.stat().st_size == EXPECTED_BYTES[(scope, split)]
        assert row["sha256"] == sha256_file(path) == EXPECTED_SHA256[(scope, split)]
        assert int(row["columns"]) == len(COLUMNS[(scope, split)])
        assert len(row["sd02_sha256"]) == 64
        assert len(row["sd03_sha256"]) == 64
        assert all(c in "0123456789abcdef" for c in row["sd02_sha256"])
        assert all(c in "0123456789abcdef" for c in row["sd03_sha256"])
        assert row["scope"] == scope and row["split"] == split
