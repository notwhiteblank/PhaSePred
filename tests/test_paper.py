from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phasepred.paper import PaperSplitError, load_paper_features, load_paper_split


def _write_workbook(path: Path, sheets: dict[str, list[str]]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, accessions in sheets.items():
            pd.DataFrame({"UniprotEntry": accessions}).to_excel(
                writer, sheet_name=sheet_name, index=False
            )


def test_load_paper_split_parses_train_and_test_sets(tmp_path: Path) -> None:
    s2 = tmp_path / "s2.xlsx"
    s3 = tmp_path / "s3.xlsx"
    _write_workbook(
        s2,
        {
            "SaPS": ["P1"],
            "PdPS": ["P2"],
            "NoPS": ["N1", "N2"],
            "hSaPS": ["HP1"],
            "hPdPS": ["HP2"],
            "hNoPS": ["HN1"],
        },
    )
    _write_workbook(
        s3,
        {
            "SaPS-test": ["P3"],
            "PdPS-test": ["P4"],
            "PS-test": ["P5"],
            "hSaPS-test": ["HP3"],
            "hPdPS-test": ["HP4"],
            "hPS-test": ["HP5"],
            "NoPS-test": ["N3"],
            "hNoPS-test": ["HN2"],
        },
    )

    split = load_paper_split(s2, s3)

    assert set(split.train["split"]) == {"train"}
    assert set(split.test["split"]) == {"test"}
    assert split.train.query("task == 'SaPS' and label == 1")["UniprotEntry"].tolist() == ["P1"]
    assert split.train.query("task == 'SaPS' and label == 0")["UniprotEntry"].tolist() == [
        "N1",
        "N2",
    ]
    ps_test = split.test.query("source_sheet == 'PS-test'")
    assert ps_test["UniprotEntry"].tolist() == ["P5", "P5"]
    assert ps_test["task"].tolist() == ["SaPS", "PdPS"]
    assert split.test.query("task == 'hPdPS' and label == 1")["UniprotEntry"].tolist() == [
        "HP4",
        "HP5",
    ]


def test_load_paper_split_requires_uniprot_entry(tmp_path: Path) -> None:
    s2 = tmp_path / "s2.xlsx"
    s3 = tmp_path / "s3.xlsx"
    with pd.ExcelWriter(s2, engine="openpyxl") as writer:
        pd.DataFrame({"wrong": ["P1"]}).to_excel(writer, sheet_name="SaPS", index=False)
    _write_workbook(s3, {"SaPS-test": ["P2"]})

    with pytest.raises(PaperSplitError):
        load_paper_split(s2, s3)


def test_load_paper_features_includes_shared_test_rows(tmp_path: Path) -> None:
    s2 = tmp_path / "s2.xlsx"
    s3 = tmp_path / "s3.xlsx"
    feature_columns = {
        "UniprotEntry": ["P1"],
        "Gene name": ["G1"],
        "Organism": ["Org"],
        "Organism ID": [1],
        "length": [4],
        "Hydropathy": [0.1],
        "FCR": [0.2],
        "IDR": [0.3],
        "LCR": [0.4],
        "PScore": [0.5],
        "PLAAC": [0.6],
        "catGRANULE": [0.7],
        "DeepCoil": [0.8],
    }
    human_columns = {
        **feature_columns,
        "Phos freq": [0.9],
        "DeepPhase": [1.0],
    }
    with pd.ExcelWriter(s2, engine="openpyxl") as writer:
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="SaPS", index=False)
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="PdPS", index=False)
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="NoPS", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hSaPS", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hPdPS", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hNoPS", index=False)
    with pd.ExcelWriter(s3, engine="openpyxl") as writer:
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="SaPS-test", index=False)
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="PdPS-test", index=False)
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="PS-test", index=False)
        pd.DataFrame(feature_columns).to_excel(writer, sheet_name="NoPS-test", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hSaPS-test", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hPdPS-test", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hPS-test", index=False)
        pd.DataFrame(human_columns).to_excel(writer, sheet_name="hNoPS-test", index=False)

    features = load_paper_features(s2, s3)

    assert "Hydropathy" in features.train.columns
    assert "Phos freq" in features.train.columns
    assert features.test.query("source_sheet == 'PS-test'")["task"].tolist() == [
        "SaPS",
        "PdPS",
    ]
    assert features.test.query("source_sheet == 'hPS-test'")["task"].tolist() == [
        "hSaPS",
        "hPdPS",
    ]
