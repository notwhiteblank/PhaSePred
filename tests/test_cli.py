from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from conftest import require_predict_requirements, require_tool
from typer.testing import CliRunner

from phasepred.cli import app

runner = CliRunner()


def _write_workbook(path: Path, sheets: dict[str, list[str]]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, accessions in sheets.items():
            pd.DataFrame({"UniprotEntry": accessions}).to_excel(
                writer, sheet_name=sheet_name, index=False
            )


def test_prepare_split_command(tmp_path: Path) -> None:
    s2 = tmp_path / "s2.xlsx"
    s3 = tmp_path / "s3.xlsx"
    _write_workbook(
        s2,
        {
            "SaPS": ["P1"],
            "PdPS": ["P2"],
            "NoPS": ["N1"],
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
            "NoPS-test": ["N2"],
            "hNoPS-test": ["HN2"],
        },
    )
    output = tmp_path / "split.csv"

    result = runner.invoke(
        app,
        [
            "prepare-split",
            "--s2",
            str(s2),
            "--s3",
            str(s3),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert set(frame["split"]) == {"train", "test"}
    assert "SaPS" in set(frame["task"])


def test_features_from_split_command(tmp_path: Path) -> None:
    split = tmp_path / "split.csv"
    sequences = tmp_path / "seqs.csv"
    output = tmp_path / "features.csv"

    pd.DataFrame(
        {
            "task": ["SaPS", "SaPS"],
            "split": ["train", "test"],
            "label": [1, 0],
            "UniprotEntry": ["P1", "P2"],
            "source_sheet": ["SaPS", "NoPS-test"],
        }
    ).to_csv(split, index=False)
    pd.DataFrame(
        {
            "UniprotEntry": ["P1", "P2"],
            "sequence": ["ACDE", "KKKK"],
        }
    ).to_csv(sequences, index=False)

    result = runner.invoke(
        app,
        [
            "features-from-split",
            "--split-table",
            str(split),
            "--sequences",
            str(sequences),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert list(frame["UniprotEntry"]) == ["P1", "P2"]
    assert "Hydropathy" in frame.columns


def test_features_from_paper_command(tmp_path: Path) -> None:
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
    output = tmp_path / "paper_features.csv"

    result = runner.invoke(
        app,
        [
            "features-from-paper",
            "--s2",
            str(s2),
            "--s3",
            str(s3),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert "Phos freq" in frame.columns
    assert frame.query("source_sheet == 'PS-test'")["task"].tolist() == ["SaPS", "PdPS"]


def test_sequences_from_ids_command_uses_cache(tmp_path: Path) -> None:
    input_path = tmp_path / "ids.csv"
    cache = tmp_path / "cache.jsonl"
    output = tmp_path / "seqs.csv"

    pd.DataFrame({"UniprotEntry": ["P12345", "P12345"]}).to_csv(input_path, index=False)
    cache.write_text(
        (
            '{"accession": "P12345", "sequence": "ACDE", "source": "fasta", '
            '"header": ">cached", "retrieved_at": null, "endpoint": null, '
            '"sequence_hash": "dummy", "length": 4}\n'
        ),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "sequences-from-ids",
            "--input",
            str(input_path),
            "--cache",
            str(cache),
            "--no-network",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert list(frame["UniprotEntry"]) == ["P12345"]
    assert list(frame["sequence"]) == ["ACDE"]


def test_train_and_predict_commands(tmp_path: Path) -> None:
    split = tmp_path / "split.csv"
    features = tmp_path / "features.csv"
    model_path = tmp_path / "model.joblib"
    output = tmp_path / "scores.csv"

    pd.DataFrame(
        {
            "task": ["SaPS", "SaPS", "SaPS", "SaPS"],
            "split": ["train", "train", "train", "train"],
            "label": [0, 0, 1, 1],
            "UniprotEntry": ["P1", "P2", "P3", "P4"],
            "source_sheet": ["NoPS", "NoPS", "SaPS", "SaPS"],
        }
    ).to_csv(split, index=False)
    pd.DataFrame(
        {
            "task": ["SaPS", "SaPS", "SaPS", "SaPS"],
            "UniprotEntry": ["P1", "P2", "P3", "P4"],
            "Hydropathy": [1.0, 1.2, -1.0, -1.2],
            "FCR": [0.1, 0.2, 0.8, 0.9],
            "IDR": [0.0, 0.1, 0.8, 0.9],
            "LCR": [0.2, 0.2, 0.7, 0.8],
            "PScore": [0.1, 0.2, 0.9, 0.8],
            "PLAAC": [0.0, 0.1, 0.8, 0.9],
            "catGRANULE": [0.1, 0.1, 0.7, 0.8],
            "DeepCoil": [0.0, 0.1, 0.8, 0.9],
        }
    ).to_csv(features, index=False)

    train_result = runner.invoke(
        app,
        [
            "train",
            "--split-table",
            str(split),
            "--features",
            str(features),
            "--task",
            "SaPS",
            "--output",
            str(model_path),
        ],
    )
    assert train_result.exit_code == 0, train_result.output

    predict_result = runner.invoke(
        app,
        [
            "predict-features",
            "--model",
            str(model_path),
            "--features",
            str(features),
            "--output",
            str(output),
        ],
    )
    assert predict_result.exit_code == 0, predict_result.output

    predictions = pd.read_csv(output)
    assert "score" in predictions.columns
    assert len(predictions) == 4
def test_features_from_fasta_command(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">sp|P12345|GENE_HUMAN example\n" + "M" + "S" * 120 + "\n", encoding="utf-8")
    output = tmp_path / "features.csv"

    result = runner.invoke(
        app,
        [
            "features-from-fasta",
            "--input",
            str(fasta),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert list(frame["UniprotEntry"]) == ["P12345"]
    assert "Hydropathy" in frame.columns
    # v2022: ESpritz-domain IDR fraction + SEG-domain LCR fraction
    assert "espritz-idr-fraction" in frame.columns
    assert "lcr-fraction" in frame.columns
    assert 0.0 <= frame.loc[0, "lcr-fraction"] <= 1.0
    idr_ok = pd.isna(frame.loc[0, "espritz-idr-fraction"]) or (
        0.0 <= frame.loc[0, "espritz-idr-fraction"] <= 1.0
    )
    assert idr_ok


def test_features_from_fasta_degrades_nan_on_missing_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing external tools must degrade espritz/lcr columns to NaN with a
    warning instead of failing (v2022 definitions)."""
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">P1\n" + "M" + "S" * 120 + "\n", encoding="utf-8")
    output = tmp_path / "features.csv"

    import phasepred.tools as tools

    def _boom(records):
        raise RuntimeError("ESpritz binary unavailable")

    monkeypatch.setattr(tools, "run_espritz", _boom)

    result = runner.invoke(
        app,
        [
            "features-from-fasta",
            "--input",
            str(fasta),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "interrupted" not in result.output.lower()
    frame = pd.read_csv(output)
    assert pd.isna(frame.loc[0, "espritz-idr-fraction"])
    # lcr still computed (SEG present)
    assert 0.0 <= frame.loc[0, "lcr-fraction"] <= 1.0


@pytest.mark.tools
def test_features_recomputed_requires_external_deepcoil(tmp_path: Path) -> None:
    """catGRANULE now comes from the `catgranule` package; the remaining
    external requirement is the --deepcoil-features table."""
    # scripts/features-recomputed computes LCR/PScore/PLAAC before it reaches
    # the --deepcoil-features check, so those three contracts must resolve.
    require_tool("seg")
    require_tool("plaac")
    require_tool("pscore")
    sequences = tmp_path / "sequences.csv"
    output = tmp_path / "features.csv"
    espritz_cache = tmp_path / "espritz.jsonl"
    pd.DataFrame({"UniprotEntry": ["P1"], "sequence": ["Q" * 90]}).to_csv(
        sequences,
        index=False,
    )
    espritz_cache.write_text(
        (
            '{"accession": "P1", "disorder_threshold": 0.5072, "espritz_idr": 0.5, '
            '"max_score": 0.8, "mean_score": 0.4, "model": "D", "residues": 90, "sw": 0}\n'
        ),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "features-recomputed",
            "--sequences",
            str(sequences),
            "--output",
            str(output),
            "--espritz-cache",
            str(espritz_cache),
        ],
    )

    assert result.exit_code != 0
    assert "--deepcoil-features is required" in result.output.replace("\n", " ")


@pytest.mark.tools
def test_features_recomputed_csv_source_requires_catgranule_features(tmp_path: Path) -> None:
    require_tool("seg")
    require_tool("plaac")
    require_tool("pscore")
    sequences = tmp_path / "sequences.csv"
    output = tmp_path / "features.csv"
    espritz_cache = tmp_path / "espritz.jsonl"
    pd.DataFrame({"UniprotEntry": ["P1"], "sequence": ["Q" * 90]}).to_csv(
        sequences,
        index=False,
    )
    espritz_cache.write_text(
        (
            '{"accession": "P1", "disorder_threshold": 0.5072, "espritz_idr": 0.5, '
            '"max_score": 0.8, "mean_score": 0.4, "model": "D", "residues": 90, "sw": 0}\n'
        ),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "features-recomputed",
            "--sequences",
            str(sequences),
            "--output",
            str(output),
            "--espritz-cache",
            str(espritz_cache),
            "--catgranule-source",
            "csv",
        ],
    )

    assert result.exit_code != 0
    assert "--catgranule-features is required" in result.output.replace("\n", " ")


def _paper_feature_row(accession: str, idr: float) -> dict[str, object]:
    return {
        "UniprotEntry": accession,
        "Gene name": accession,
        "Organism": "Org",
        "Organism ID": 1,
        "length": 5,
        "Hydropathy": idr,
        "FCR": idr,
        "IDR": idr,
        "LCR": idr,
        "PScore": idr,
        "PLAAC": idr,
        "catGRANULE": idr,
        "DeepCoil": idr,
    }


@pytest.mark.tools
def test_predict_no_product_flag_uses_single_product(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Product B is gone: predict has no --product and scores with the
    default Product A artifact (or --models-dir)."""
    require_predict_requirements()
    from phasepred.features import BASE_FEATURE_COLUMNS
    from phasepred.models import train_model_artifact

    model_root = tmp_path / "models"
    mode_dir = model_root / "SaPS"
    mode_dir.mkdir(parents=True)

    feature_rows = [
        {"UniprotEntry": "P1", "label": 1, "length": 5.0, **{c: 0.5 for c in BASE_FEATURE_COLUMNS}},
        {"UniprotEntry": "P2", "label": 0, "length": 5.0, **{c: 0.1 for c in BASE_FEATURE_COLUMNS}},
    ]
    fake_features = pd.DataFrame(feature_rows)
    import joblib

    for i in range(2):
        artifact = train_model_artifact(fake_features, label_column="label", task="SaPS")
        joblib.dump(artifact, mode_dir / f"8f_model_{i}.joblib")

    monkeypatch.setattr(
        "phasepred.predictor.compute_all_features",
        lambda records, include_human=False: fake_features.drop(columns=["label"]),
    )

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\nACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDE\n", encoding="utf-8")
    output = tmp_path / "out.csv"

    result = runner.invoke(
        app,
        [
            "predict",
            "--fasta",
            str(fasta),
            "--mode",
            "SaPS",
            "--models-dir",
            str(model_root),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    frame = pd.read_csv(output)
    assert "stale" in frame.columns
    assert frame["stale"].all()
    assert "STALE" in result.output


@pytest.mark.tools
def test_predict_stale_warning_suppressed_fresh_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    require_predict_requirements()
    import json

    import joblib

    from phasepred.features import BASE_FEATURE_COLUMNS
    from phasepred.models import FEATURE_DEFINITIONS_VERSION, train_model_artifact

    model_root = tmp_path / "models"
    mode_dir = model_root / "SaPS"
    mode_dir.mkdir(parents=True)
    mode_dir.joinpath("manifest.json").write_text(
        json.dumps({"feature_definitions_version": FEATURE_DEFINITIONS_VERSION}),
        encoding="utf-8",
    )

    row1 = {
        "UniprotEntry": "P1", "label": 1, "length": 48.0, **{c: 0.5 for c in BASE_FEATURE_COLUMNS}
    }
    row2 = {
        "UniprotEntry": "P2", "label": 0, "length": 48.0, **{c: 0.1 for c in BASE_FEATURE_COLUMNS}
    }
    fake_features = pd.DataFrame([row1, row2])
    for i in range(2):
        joblib.dump(
            train_model_artifact(fake_features, label_column="label", task="SaPS"),
            mode_dir / f"8f_model_{i}.joblib",
        )

    monkeypatch.setattr(
        "phasepred.predictor.compute_all_features",
        lambda records, include_human=False: fake_features.drop(columns=["label"]),
    )

    fasta = tmp_path / "in.fasta"
    fasta.write_text(
        ">P1\nACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDEACDE\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.csv"

    result = runner.invoke(
        app,
        [
            "predict",
            "--fasta",
            str(fasta),
            "--mode",
            "SaPS",
            "--models-dir",
            str(model_root),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "STALE" not in result.output
    frame = pd.read_csv(output)
    assert not frame["stale"].any()
