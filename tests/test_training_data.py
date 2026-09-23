"""Format contract for data/processed/*.tsv (E1).

These tests use a miniature stand-in for ``load_paper_features()`` output, so
they need neither the private supplementary workbooks nor any external tool.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS
from phasepred.training_data import (
    COLUMNS,
    DATA_ACCESSION_DATE,
    EXPECTED_BYTES,
    EXPECTED_ROWS,
    EXPECTED_SHA256,
    FILENAMES,
    FLOAT_COLUMNS,
    INT_COLUMNS,
    MANIFEST_COLUMNS,
    SCOPES,
    SCRIPT_VERSION,
    SPLITS,
    TASKS_BY_SCOPE,
    TEXT_COLUMNS,
    TrainingDataSchemaError,
    dtype_map,
    load_task_frame,
    manifest_path,
    manifest_row,
    processed_dir,
    read_manifest,
    read_table,
    read_tsv,
    select_and_cast,
    sha256_file,
    table_path,
    validate_expected_size,
    validate_schema,
    write_manifest,
    write_tsv,
)

ACC = [f"P0000{i}" for i in range(1, 9)]


def source_frame() -> pd.DataFrame:
    """Miniature stand-in for ``load_paper_features().train``.

    Mirrors the real sheets' quirks: unnamed extra columns (``NoPS`` carries 20
    columns of which 13 are named), missing values in text/int/float columns,
    empty strings that must not be confused with missing values, and rows
    ordered by task then sheet.
    """
    return pd.DataFrame(
        {
            "UniprotEntry": ACC,
            "task": ["SaPS", "SaPS", "SaPS", "PdPS", "PdPS", "hSaPS", "hSaPS", "hPdPS"],
            "split": ["train"] * 8,
            "label": [1, 0, 0, 1, 0, 1, 0, 1],
            "source_sheet": [
                "SaPS", "NoPS", "NoPS", "PdPS", "NoPS", "hSaPS", "hNoPS", "hPdPS",
            ],
            "Gene name": ["A1", "A2", None, "B1", "B2", "C1", "", "C3"],
            "Organism": [
                "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Homo sapiens",
                "Mus musculus", "Homo sapiens", "Homo sapiens", "Homo sapiens",
            ],
            "Organism ID": [9606, 10090, None, 9606, 10090, 9606, 9606, 9606],
            "length": [1270, 300, 50, 400, 12, 999, 64, 526],
            "Hydropathy": [0.4623534558, -0.1, 0.0, 0.25, np.nan, 0.5, 0.125, -0.25],
            "FCR": [0.2251968504, 0.1, 0.0, 0.3, 0.4, 0.5, 0.6, 0.7],
            "IDR": [0.0, 0.5, 1e-05, 0.25, 0.75, 0.125, 0.0, 1.0],
            "LCR": [0.1125984252, 0.0, 0.2, 0.1, 0.0, 0.3, np.nan, 0.4],
            "PScore": [8.96, 1.0, 0.0, 2.5, 3.5, 4.5, 5.5, 6.5],
            "PLAAC": [0.224, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            "catGRANULE": [1.43649, -0.5, 0.0, 1.0, 2.0, 3.0, -1.0, 5.75],
            "DeepCoil": [0, 1, 0, 0, 1, 1, 0, 1],
            "Phos freq": [0.0440944882, np.nan, 0.0, 0.1, 0.2, 0.3, np.nan, 0.4],
            "DeepPhase": [0.9989901, 0.5, np.nan, 0.25, 0.75, 0.125, 0.0, 0.99],
            "Unnamed: 13": ["junk"] * 8,
            "Unnamed: 14": [None] * 8,
        }
    )


def source_test_frame() -> pd.DataFrame:
    """Stand-in for ``load_paper_features().test``: adds the leading Source."""
    frame = source_frame().head(3).copy()
    frame["task"] = ["SaPS", "SaPS", "PdPS"]
    frame["split"] = "test"
    frame["source_sheet"] = ["SaPS-test", "PS-test", "NoPS-test"]
    frame["Source"] = ["PhaSepDB", "PhaSepDB", None]
    return frame


def assert_text_equal(got: pd.Series, want: pd.Series, column: str) -> None:
    """Compare two string-dtype series without tripping on ``pd.NA == pd.NA``."""
    assert got.isna().tolist() == want.isna().tolist(), f"{column}: NA mask differs"
    assert got.dropna().tolist() == want.dropna().tolist(), f"{column}: text differs"


def assert_numeric_equal(got: pd.Series, want: pd.Series, column: str) -> None:
    gv = pd.to_numeric(got, errors="coerce").to_numpy(dtype="float64")
    wv = pd.to_numeric(want, errors="coerce").to_numpy(dtype="float64")
    assert np.array_equal(gv, wv, equal_nan=True), f"{column}: values differ"


# --------------------------------------------------------------- constants


def test_constant_maps_share_the_same_keys():
    keys = {(scope, split) for scope in SCOPES for split in SPLITS}
    assert set(COLUMNS) == keys
    assert set(FILENAMES) == keys
    assert set(EXPECTED_ROWS) == keys
    assert set(EXPECTED_BYTES) == keys
    assert set(EXPECTED_SHA256) <= keys


def test_column_counts_match_the_audited_widths():
    assert len(COLUMNS[("base", "train")]) == 17
    assert len(COLUMNS[("base", "test")]) == 18
    assert len(COLUMNS[("human", "train")]) == 19
    assert len(COLUMNS[("human", "test")]) == 20


def test_test_tables_prepend_the_source_column():
    assert COLUMNS[("base", "test")][0] == "Source"
    assert COLUMNS[("human", "test")][0] == "Source"
    assert COLUMNS[("base", "test")][1:] == COLUMNS[("base", "train")]
    assert COLUMNS[("human", "test")][1:] == COLUMNS[("human", "train")]


def test_feature_columns_come_from_the_v2022_definitions():
    base = COLUMNS[("base", "train")]
    human = COLUMNS[("human", "train")]
    assert base[-len(BASE_FEATURE_COLUMNS):] == list(BASE_FEATURE_COLUMNS)
    assert human[-len(HUMAN_FEATURE_COLUMNS):] == list(HUMAN_FEATURE_COLUMNS)


def test_filenames_encode_scope_and_split():
    assert FILENAMES[("base", "train")] == "chen2022_s2s3_base-features_train.tsv"
    assert FILENAMES[("base", "test")] == "chen2022_s2s3_base-features_test.tsv"
    assert FILENAMES[("human", "train")] == "chen2022_s2s3_human-features_train.tsv"
    assert FILENAMES[("human", "test")] == "chen2022_s2s3_human-features_test.tsv"


def test_expected_sizes_are_the_audited_values():
    assert EXPECTED_ROWS[("base", "train")] == 96658
    assert EXPECTED_ROWS[("base", "test")] == 24416
    assert EXPECTED_ROWS[("human", "train")] == 17757
    assert EXPECTED_ROWS[("human", "test")] == 4540
    assert sum(EXPECTED_ROWS.values()) == 143371
    assert EXPECTED_BYTES[("base", "train")] == 13399739
    assert EXPECTED_BYTES[("base", "test")] == 3508476
    assert EXPECTED_BYTES[("human", "train")] == 2800875
    assert EXPECTED_BYTES[("human", "test")] == 736518


def test_tasks_by_scope_partitions_the_four_modes():
    assert TASKS_BY_SCOPE["base"] == ("SaPS", "PdPS")
    assert TASKS_BY_SCOPE["human"] == ("hSaPS", "hPdPS")


def test_column_kind_sets_are_disjoint_and_cover_the_schema():
    assert not (INT_COLUMNS & FLOAT_COLUMNS)
    assert not (INT_COLUMNS & TEXT_COLUMNS)
    assert not (FLOAT_COLUMNS & TEXT_COLUMNS)
    for columns in COLUMNS.values():
        for column in columns:
            assert column in INT_COLUMNS | FLOAT_COLUMNS | TEXT_COLUMNS, column


def test_dtype_map_is_complete_and_stable():
    mapping = dtype_map("base", "train")
    assert set(mapping) == set(COLUMNS[("base", "train")])
    assert mapping["length"] == "Int64"
    assert mapping["DeepCoil"] == "Int64"
    assert mapping["label"] == "Int64"
    assert mapping["Organism ID"] == "Int64"
    assert mapping["Hydropathy"] == "float64"
    assert mapping["UniprotEntry"] == "string"
    assert dtype_map("human", "test")["Source"] == "string"


# ------------------------------------------------------------------- paths


def test_processed_dir_defaults_under_the_data_root():
    assert processed_dir().parts[-2:] == ("data", "processed")


def test_processed_dir_root_override_is_used_verbatim(tmp_path):
    assert processed_dir(tmp_path) == tmp_path
    name = FILENAMES[("base", "train")]
    assert table_path("base", "train", root=tmp_path) == tmp_path / name
    assert manifest_path(root=tmp_path) == tmp_path / "MANIFEST.tsv"


def test_unknown_scope_or_split_is_rejected(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="scope"):
        table_path("nope", "train", root=tmp_path)
    with pytest.raises(TrainingDataSchemaError, match="split"):
        table_path("base", "nope", root=tmp_path)


# --------------------------------------------------------- select_and_cast


def test_select_and_cast_drops_unnamed_extras_and_fixes_order():
    out = select_and_cast(source_frame(), "base", "train")
    assert list(out.columns) == COLUMNS[("base", "train")]
    assert "Unnamed: 13" not in out.columns
    assert "Unnamed: 14" not in out.columns
    assert "Phos freq" not in out.columns
    assert "DeepPhase" not in out.columns


def test_select_and_cast_keeps_human_features_only_for_human_scope():
    out = select_and_cast(source_frame(), "human", "train")
    assert list(out.columns) == COLUMNS[("human", "train")]
    assert "Phos freq" in out.columns
    assert "DeepPhase" in out.columns


def test_select_and_cast_prepends_source_for_test_tables():
    out = select_and_cast(source_test_frame(), "base", "test")
    assert list(out.columns) == COLUMNS[("base", "test")]
    assert out.columns[0] == "Source"
    assert out["Source"].isna().tolist() == [False, False, True]
    assert out["Source"].dropna().tolist() == ["PhaSepDB", "PhaSepDB"]


def test_select_and_cast_filters_scope_and_split_without_reordering():
    base = select_and_cast(source_frame(), "base", "train")
    assert base["UniprotEntry"].tolist() == ACC[:5]
    human = select_and_cast(source_frame(), "human", "train")
    assert human["UniprotEntry"].tolist() == ACC[5:]
    assert select_and_cast(source_frame(), "base", "test").empty


def test_select_and_cast_applies_the_declared_dtypes():
    out = select_and_cast(source_frame(), "human", "train")
    for column, expected in (("length", "Int64"), ("DeepCoil", "Int64"),
                             ("Organism ID", "Int64"), ("Hydropathy", "float64"),
                             ("Gene name", "string")):
        assert str(out[column].dtype) == expected, column


def test_select_and_cast_normalizes_empty_strings_to_missing():
    out = select_and_cast(source_frame(), "human", "train")
    # human scope keeps P00006 ("C1"), P00007 ("") and P00008 ("C3").
    assert out["Gene name"].isna().tolist() == [False, True, False]


def test_select_and_cast_rejects_a_frame_missing_schema_columns():
    frame = source_frame().drop(columns=["LCR"])
    with pytest.raises(TrainingDataSchemaError, match="LCR"):
        select_and_cast(frame, "base", "train")


# ------------------------------------------------------------- round trip


def test_tsv_round_trip_is_value_identical(tmp_path):
    want = select_and_cast(source_frame(), "human", "train")
    path = table_path("human", "train", root=tmp_path)
    write_tsv(want, path)
    got = read_tsv(path, "human", "train")
    assert list(got.columns) == list(want.columns)
    assert len(got) == len(want)
    for column in want.columns:
        if column in INT_COLUMNS or column in FLOAT_COLUMNS:
            assert_numeric_equal(got[column], want[column], column)
        else:
            assert_text_equal(got[column], want[column], column)


def test_float_round_trip_is_bit_exact(tmp_path):
    want = select_and_cast(source_frame(), "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    got = read_tsv(path, "base", "train")
    for column in sorted(FLOAT_COLUMNS & set(want.columns)):
        want_bytes = want[column].to_numpy(dtype="float64").tobytes()
        got_bytes = got[column].to_numpy(dtype="float64").tobytes()
        assert want_bytes == got_bytes, column


def test_float_round_trip_survives_seventeen_digit_values(tmp_path):
    frame = source_frame()
    frame.loc[0, "Hydropathy"] = 0.044094488200000004
    frame.loc[1, "Hydropathy"] = 0.46235345580000003
    want = select_and_cast(frame, "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    got = read_tsv(path, "base", "train")
    want_bytes = want["Hydropathy"].to_numpy(dtype="float64").tobytes()
    got_bytes = got["Hydropathy"].to_numpy(dtype="float64").tobytes()
    assert got_bytes == want_bytes


def test_integers_are_written_without_a_decimal_point(tmp_path):
    want = select_and_cast(source_frame(), "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    text = path.read_text(encoding="utf-8")
    assert "\t1270\t" in text
    assert "1270.0" not in text


def test_missing_values_are_written_as_empty_fields(tmp_path):
    want = select_and_cast(source_frame(), "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    header = lines[0].split("\t")
    row_p3 = dict(zip(header, lines[3].split("\t"), strict=True))
    assert row_p3["Gene name"] == ""
    assert row_p3["Organism ID"] == ""
    row_p5 = dict(zip(header, lines[5].split("\t"), strict=True))
    assert row_p5["Hydropathy"] == ""
    assert "nan" not in lines[5].lower()
    for sentinel in ("NaN", "<NA>", "NULL", "None"):
        assert sentinel not in text, sentinel


def test_text_that_looks_like_a_na_sentinel_survives(tmp_path):
    frame = source_frame()
    frame.loc[0, "Gene name"] = "NA"
    frame.loc[1, "Organism"] = "N/A"
    frame.loc[2, "Organism"] = "null"
    want = select_and_cast(frame, "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    got = read_tsv(path, "base", "train")
    assert got["Gene name"].iloc[0] == "NA"
    assert got["Organism"].iloc[1] == "N/A"
    assert got["Organism"].iloc[2] == "null"


def test_output_uses_unix_newlines_and_utf8(tmp_path):
    want = select_and_cast(source_frame(), "base", "train")
    path = table_path("base", "train", root=tmp_path)
    write_tsv(want, path)
    raw = path.read_bytes()
    assert b"\r\n" not in raw
    assert raw.endswith(b"\n")
    raw.decode("utf-8")


def test_write_tsv_creates_missing_parents(tmp_path):
    path = tmp_path / "nested" / "deeper" / FILENAMES[("base", "train")]
    write_tsv(select_and_cast(source_frame(), "base", "train"), path)
    assert path.is_file()


# -------------------------------------------------------------- validation


def test_validate_schema_accepts_a_well_formed_table():
    validate_schema(select_and_cast(source_frame(), "base", "train"), "base", "train")


def test_validate_schema_rejects_reordered_columns():
    frame = select_and_cast(source_frame(), "base", "train")
    with pytest.raises(TrainingDataSchemaError, match="column mismatch"):
        validate_schema(frame[list(reversed(frame.columns))], "base", "train")


def test_validate_schema_rejects_an_out_of_scope_task():
    frame = select_and_cast(source_frame(), "base", "train")
    frame.loc[0, "task"] = "hSaPS"
    with pytest.raises(TrainingDataSchemaError, match="out-of-scope"):
        validate_schema(frame, "base", "train")


def test_validate_schema_rejects_a_mixed_split_column():
    frame = select_and_cast(source_frame(), "base", "train")
    frame.loc[0, "split"] = "test"
    with pytest.raises(TrainingDataSchemaError, match="split"):
        validate_schema(frame, "base", "train")


def test_validate_schema_rejects_a_bad_label():
    frame = select_and_cast(source_frame(), "base", "train")
    frame.loc[0, "label"] = 2
    with pytest.raises(TrainingDataSchemaError, match="label"):
        validate_schema(frame, "base", "train")


def test_validate_schema_rejects_duplicates_within_a_task_and_split():
    frame = select_and_cast(source_frame(), "base", "train")
    duplicated = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    with pytest.raises(TrainingDataSchemaError, match="duplicated"):
        validate_schema(duplicated, "base", "train")


def test_read_tsv_rejects_a_broken_header(tmp_path):
    path = table_path("base", "train", root=tmp_path)
    write_tsv(select_and_cast(source_frame(), "base", "train"), path)
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    header[header.index("LCR")] = "LCR_renamed"
    lines[0] = "\t".join(header)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(TrainingDataSchemaError, match="column mismatch"):
        read_tsv(path, "base", "train")


def test_read_table_reports_a_missing_file_actionably(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="build_training_tsv"):
        read_table("base", "train", root=tmp_path)


def test_load_task_frame_filters_and_preserves_order(tmp_path):
    frame = select_and_cast(source_frame(), "base", "train")
    write_tsv(frame, table_path("base", "train", root=tmp_path))
    saps = load_task_frame("base", "train", "SaPS", root=tmp_path)
    assert saps["UniprotEntry"].tolist() == ACC[:3]
    assert saps.index.tolist() == [0, 1, 2]
    pdps = load_task_frame("base", "train", "PdPS", root=tmp_path)
    assert pdps["UniprotEntry"].tolist() == ACC[3:5]
    assert pdps.index.tolist() == [0, 1]


def test_load_task_frame_rejects_a_task_outside_the_scope(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="does not belong"):
        load_task_frame("base", "train", "hSaPS", root=tmp_path)


def test_load_task_frame_rejects_an_unknown_scope(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="scope"):
        load_task_frame("nope", "train", "SaPS", root=tmp_path)


# ---------------------------------------------------------------- manifest


def test_manifest_row_has_no_wall_clock_field(tmp_path):
    path = table_path("base", "train", root=tmp_path)
    write_tsv(select_and_cast(source_frame(), "base", "train"), path)
    row = manifest_row(path, "base", "train", rows=5,
                       sd02_sha256="aa" * 32, sd03_sha256="bb" * 32)
    assert set(row) == set(MANIFEST_COLUMNS)
    assert row["file"] == FILENAMES[("base", "train")]
    assert row["scope"] == "base"
    assert row["split"] == "train"
    assert row["rows"] == 5
    assert row["columns"] == 17
    assert row["bytes"] == path.stat().st_size
    assert row["sha256"] == sha256_file(path)
    assert row["script_version"] == SCRIPT_VERSION
    assert row["data_accession_date"] == DATA_ACCESSION_DATE


def test_manifest_write_is_deterministic(tmp_path):
    path = table_path("base", "train", root=tmp_path)
    write_tsv(select_and_cast(source_frame(), "base", "train"), path)
    row = manifest_row(path, "base", "train", rows=5,
                       sd02_sha256="aa" * 32, sd03_sha256="bb" * 32)
    first, second = tmp_path / "m1.tsv", tmp_path / "m2.tsv"
    write_manifest([row], first)
    write_manifest([row], second)
    assert first.read_bytes() == second.read_bytes()
    assert sha256_file(first) == sha256_file(second)


def test_manifest_round_trips_through_read_manifest(tmp_path):
    path = table_path("base", "train", root=tmp_path)
    write_tsv(select_and_cast(source_frame(), "base", "train"), path)
    row = manifest_row(path, "base", "train", rows=5,
                       sd02_sha256="aa" * 32, sd03_sha256="bb" * 32)
    write_manifest([row], manifest_path(root=tmp_path))
    got = read_manifest(root=tmp_path)
    assert list(got.columns) == MANIFEST_COLUMNS
    assert got.iloc[0]["file"] == FILENAMES[("base", "train")]
    assert got.iloc[0]["sha256"] == row["sha256"]


def test_read_manifest_reports_a_missing_file(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="Manifest missing"):
        read_manifest(root=tmp_path)


def test_validate_expected_size_flags_tampering(tmp_path, monkeypatch):
    scope, split = "base", "train"
    key = (scope, split)
    path = table_path(scope, split, root=tmp_path)
    write_tsv(select_and_cast(source_frame(), scope, split), path)
    monkeypatch.setitem(EXPECTED_ROWS, key, 5)
    monkeypatch.setitem(EXPECTED_BYTES, key, path.stat().st_size)
    monkeypatch.setitem(EXPECTED_SHA256, key, sha256_file(path))
    validate_expected_size(scope, split, root=tmp_path)
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(TrainingDataSchemaError, match="bytes"):
        validate_expected_size(scope, split, root=tmp_path)


def test_validate_expected_size_flags_a_row_count_drift(tmp_path, monkeypatch):
    scope, split = "base", "train"
    key = (scope, split)
    path = table_path(scope, split, root=tmp_path)
    write_tsv(select_and_cast(source_frame(), scope, split), path)
    monkeypatch.setitem(EXPECTED_ROWS, key, 6)
    monkeypatch.setitem(EXPECTED_BYTES, key, path.stat().st_size)
    monkeypatch.setitem(EXPECTED_SHA256, key, sha256_file(path))
    with pytest.raises(TrainingDataSchemaError, match="rows"):
        validate_expected_size(scope, split, root=tmp_path)


def test_validate_expected_size_flags_a_digest_drift(tmp_path, monkeypatch):
    scope, split = "base", "train"
    key = (scope, split)
    path = table_path(scope, split, root=tmp_path)
    write_tsv(select_and_cast(source_frame(), scope, split), path)
    monkeypatch.setitem(EXPECTED_ROWS, key, 5)
    monkeypatch.setitem(EXPECTED_BYTES, key, path.stat().st_size)
    monkeypatch.setitem(EXPECTED_SHA256, key, "cc" * 32)
    with pytest.raises(TrainingDataSchemaError, match="sha256"):
        validate_expected_size(scope, split, root=tmp_path)


def test_validate_expected_size_reports_a_missing_table(tmp_path):
    with pytest.raises(TrainingDataSchemaError, match="missing"):
        validate_expected_size("base", "train", root=tmp_path)


def test_sha256_file_matches_hashlib(tmp_path):
    import hashlib

    path = tmp_path / "blob.bin"
    path.write_bytes(b"phasepred" * 1000)
    assert sha256_file(path) == hashlib.sha256(path.read_bytes()).hexdigest()
