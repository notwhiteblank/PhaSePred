"""Unit tests for the three-layer validation module (E3, Task 1).

Every case here uses synthetic metrics, synthetic training tables written to
``tmp_path`` and synthetic model directories.  The module therefore runs in a
public checkout: no private supplementary workbook, no gitignored sequence
cache and no external tool is required.

Each case documents the reverse proof that makes it non-vacuous: the property
it asserts is one a plausible implementation could get wrong, and the paired
case exercises that failure.  The train/test accession overlap test is the
canonical example: the positive half flags a planted accession and the negative
half removes it, so a check that always passes or always fails turns one of the
two red.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xgboost
from sklearn.metrics import roc_auc_score

from phasepred.models import FEATURE_DEFINITIONS_VERSION
from phasepred.training import (
    N_FOLDS,
    N_MODELS,
    NEG_RATIO,
    SEED,
    TASKS,
    feature_columns,
)
from phasepred.training_data import COLUMNS, table_path, write_tsv
from phasepred.validation import (
    PAPER_TARGETS,
    PAPER_TOLERANCE,
    CheckResult,
    ValidationError,
    _accessions,
    _auc,
    _relative,
    auc_layer,
    consistency_layer,
    exit_code_for,
    format_report_table,
    leakage_layer,
    recomputed_auc_layer,
    results_to_report,
)

# --------------------------------------------------------------- fixtures


def _tables() -> dict[tuple[str, str], list[dict[str, str]]]:
    return {
        ("base", "train"): [
            {"UniprotEntry": "SA_T1", "task": "SaPS", "split": "train",
             "source_sheet": "SaPS", "label": "1"},
            {"UniprotEntry": "PD_T1", "task": "PdPS", "split": "train",
             "source_sheet": "PdPS", "label": "1"},
            {"UniprotEntry": "NEG_B1", "task": "SaPS", "split": "train",
             "source_sheet": "NoPS", "label": "0"},
            {"UniprotEntry": "NEG_B2", "task": "PdPS", "split": "train",
             "source_sheet": "NoPS", "label": "0"},
        ],
        ("base", "test"): [
            {"UniprotEntry": "SA_E1", "task": "SaPS", "split": "test",
             "source_sheet": "SaPS-test", "label": "1"},
            {"UniprotEntry": "PD_E1", "task": "PdPS", "split": "test",
             "source_sheet": "PdPS-test", "label": "1"},
            {"UniprotEntry": "NEG_BE1", "task": "SaPS", "split": "test",
             "source_sheet": "NoPS-test", "label": "0"},
            {"UniprotEntry": "NEG_BE2", "task": "PdPS", "split": "test",
             "source_sheet": "NoPS-test", "label": "0"},
        ],
        ("human", "train"): [
            {"UniprotEntry": "HSA_T1", "task": "hSaPS", "split": "train",
             "source_sheet": "hSaPS", "label": "1"},
            {"UniprotEntry": "HPD_T1", "task": "hPdPS", "split": "train",
             "source_sheet": "hPdPS", "label": "1"},
            {"UniprotEntry": "NEG_H1", "task": "hSaPS", "split": "train",
             "source_sheet": "hNoPS", "label": "0"},
            {"UniprotEntry": "NEG_H2", "task": "hPdPS", "split": "train",
             "source_sheet": "hNoPS", "label": "0"},
        ],
        ("human", "test"): [
            {"UniprotEntry": "HSA_E1", "task": "hSaPS", "split": "test",
             "source_sheet": "hSaPS-test", "label": "1"},
            {"UniprotEntry": "HPD_E1", "task": "hPdPS", "split": "test",
             "source_sheet": "hPdPS-test", "label": "1"},
            {"UniprotEntry": "NEG_HE1", "task": "hSaPS", "split": "test",
             "source_sheet": "hNoPS-test", "label": "0"},
            {"UniprotEntry": "NEG_HE2", "task": "hPdPS", "split": "test",
             "source_sheet": "hNoPS-test", "label": "0"},
        ],
    }


def _write_all(root: Path, tables: dict[tuple[str, str], list[dict[str, str]]]) -> None:
    for (scope, split), rows in tables.items():
        columns = COLUMNS[(scope, split)]
        frame = pd.DataFrame(
            [{**{column: "" for column in columns}, **row} for row in rows],
            columns=columns,
        )
        write_tsv(frame, table_path(scope, split, root=root))


def _by_name(results: list[CheckResult], name: str) -> CheckResult:
    matches = [result for result in results if result.name == name]
    assert len(matches) == 1, f"expected exactly one {name!r}, got {[r.name for r in results]}"
    return matches[0]


def _manifest(task: str, **overrides: object) -> dict[str, object]:
    data: dict[str, object] = {
        "task": task,
        "feature_columns": feature_columns(task),
        "feature_definitions_version": FEATURE_DEFINITIONS_VERSION,
        "protocol": "paper-default",
        "n_models": N_MODELS,
        "n_neg_ratio": NEG_RATIO,
        "n_folds": N_FOLDS,
        "seed": SEED,
        "xgboost_version": xgboost.__version__,
    }
    data.update(overrides)
    return data


def _make_artifact(
    models_dir: Path, task: str, *, n_files: int = N_MODELS, **overrides: object
) -> None:
    directory = models_dir / task
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    (directory / "manifest.json").write_text(
        json.dumps(_manifest(task, **overrides)), encoding="utf-8"
    )
    for index in range(n_files):
        (directory / f"8f_model_{index}.joblib").touch()


def _make_all(models_dir: Path) -> None:
    for task in TASKS:
        _make_artifact(models_dir, task)


def _result(
    layer: str = "paper",
    name: str = "row",
    verdict: str = "fail",
    hard: bool = False,
    measured: object = 0.0,
    threshold: object = 0.0,
    attribution: str = "",
) -> CheckResult:
    return CheckResult(
        layer=layer,
        name=name,
        verdict=verdict,
        hard=hard,
        measured=measured,
        threshold=threshold,
        attribution=attribution,
    )


def _walk_keys(value: object) -> set[str]:
    keys: set[str] = set()
    if isinstance(value, dict):
        for key, item in value.items():
            keys.add(key)
            keys |= _walk_keys(item)
    elif isinstance(value, list):
        for item in value:
            keys |= _walk_keys(item)
    return keys


# -------------------------------------------------------------- AUC layer

def test_auc_layer_passes_a_row_within_tolerance():
    row = _by_name(auc_layer({"SaPS": {"cv_auc": 0.862 - 0.005}}), "SaPS-8")
    assert row.verdict == "pass"
    assert row.hard is False


def test_auc_layer_fails_a_row_beyond_tolerance():
    row = _by_name(auc_layer({"hPdPS": {"cv_auc": 0.827 - 0.0108}}), "hPdPS-10")
    assert row.verdict == "fail"
    assert row.hard is False
    assert row.attribution
    assert "S5" in row.attribution
    assert "SEED=42" in row.attribution


def test_auc_layer_uses_per_model_not_ensemble():
    metrics = {"hPdPS": {"cv_auc": 0.827 - 0.05, "cv_auc_ensemble": 0.827 + 0.001}}
    row = _by_name(auc_layer(metrics), "hPdPS-10")
    assert row.verdict == "fail"
    assert row.measured == pytest.approx(0.827 - 0.05)


def test_auc_layer_covers_the_four_paper_rows():
    assert set(PAPER_TARGETS) == {
        ("SaPS", 8),
        ("PdPS", 8),
        ("hSaPS", 10),
        ("hPdPS", 10),
    }
    assert PAPER_TOLERANCE == 0.01


def test_auc_layer_skips_a_task_absent_from_metrics():
    row = _by_name(auc_layer({"SaPS": {"cv_auc": 0.862}}), "hPdPS-10")
    assert row.verdict == "skip"
    assert row.attribution


def test_auc_helper_agrees_with_sklearn_on_a_hand_computed_case():
    labels = np.array([0, 0, 1, 1])
    scores = np.array([0.1, 0.4, 0.35, 0.8])
    # Positive/negative pairs correctly ordered: (0.35>0.1), (0.8>0.1), (0.8>0.4)
    # out of four pairs, so the hand value is 0.75.
    assert _auc(labels, scores) == pytest.approx(0.75)
    assert _auc(labels, scores) == pytest.approx(roc_auc_score(labels, scores))


def test_auc_helper_returns_nan_for_a_single_class():
    value = _auc(np.array([1, 1, 1]), np.array([0.1, 0.5, 0.9]))
    assert np.isnan(value)


def test_recomputed_auc_layer_requires_a_metrics_file(tmp_path):
    with pytest.raises(ValidationError):
        recomputed_auc_layer(metrics_path=tmp_path / "missing_metrics.json")


def test_recomputed_auc_layer_sends_the_frozen_progress_to_stderr(
    tmp_path, monkeypatch, capsys
):
    """Ruling 8: the progress print of the frozen body must not reach stdout.

    ``run_cross_validation`` prints one line per fold from inside its frozen
    protocol body (``training.py:158``), so the only place to redirect it is the
    call site in ``recomputed_auc_layer``. The fake below prints exactly what the
    frozen body prints; removing the ``redirect_stdout`` wrapper puts the line on
    stdout and fails this test. ``ensemble_score`` and ``joblib.load`` are
    stubbed so the module's own recomputation loop is still the code under test.
    """
    import phasepred.validation as validation

    _write_all(tmp_path, _tables())
    models_dir = tmp_path / "models"
    _make_all(models_dir)

    def noisy_cross_validation(*args, **kwargs):
        print("    fold 1/5: PROGRESS noise that used to pollute stdout")
        return 0.86, [0.86], 0.86, [0.86]

    monkeypatch.setattr(validation, "run_cross_validation", noisy_cross_validation)
    monkeypatch.setattr(validation.joblib, "load", lambda path: object())
    monkeypatch.setattr(validation, "ensemble_score", lambda models, X: np.zeros(len(X)))

    metrics = {
        task: {"cv_auc": 0.86, "cv_auc_ensemble": 0.86, "test_auc_S3": 0.5}
        for task in TASKS
    }
    metrics_path = tmp_path / "metrics.json"
    metrics_path.write_text(json.dumps(metrics), encoding="utf-8")

    recomputed_auc_layer(root=tmp_path, models_dir=models_dir, metrics_path=metrics_path)

    captured = capsys.readouterr()
    assert "PROGRESS noise" in captured.err
    assert captured.out == ""


def test_accessions_strips_surrounding_whitespace():
    frame = pd.DataFrame({"UniprotEntry": ["  P12345  ", "P12345", " Q9Y6  "]})
    assert _accessions(frame) == {"P12345", "Q9Y6"}


# ---------------------------------------------------------- leakage layer


@pytest.mark.parametrize("overlap", [True, False])
def test_leakage_flags_a_train_test_accession_overlap(tmp_path, overlap):
    tables = _tables()
    if overlap:
        tables[("base", "train")].append(
            {"UniprotEntry": "LEAK1", "task": "SaPS", "split": "train",
             "source_sheet": "SaPS", "label": "1"}
        )
        tables[("base", "test")].append(
            {"UniprotEntry": "LEAK1", "task": "SaPS", "split": "test",
             "source_sheet": "SaPS-test", "label": "1"}
        )
    _write_all(tmp_path, tables)
    row = _by_name(leakage_layer(root=tmp_path), "SaPS: train/test accession overlap")
    assert row.hard is True
    if overlap:
        assert row.verdict == "fail"
        assert "LEAK1" in row.measured["overlap"]
        assert exit_code_for([row]) == 1
    else:
        assert row.verdict == "pass"
        assert exit_code_for([row]) == 0


@pytest.mark.parametrize("label", ["1", "0"])
def test_leakage_flags_a_positive_across_the_split(tmp_path, label):
    tables = _tables()
    tables[("base", "train")].append(
        {"UniprotEntry": "XAC", "task": "SaPS", "split": "train",
         "source_sheet": "SaPS", "label": label}
    )
    tables[("base", "test")].append(
        {"UniprotEntry": "XAC", "task": "SaPS", "split": "test",
         "source_sheet": "SaPS-test", "label": label}
    )
    _write_all(tmp_path, tables)
    results = leakage_layer(root=tmp_path)
    row = _by_name(results, "SaPS: positive in both train and test")
    assert row.hard is True
    if label == "1":
        assert row.verdict == "fail"
        assert "XAC" in row.measured["overlap"]
        assert exit_code_for([row]) == 1
    else:
        assert row.verdict == "pass"
    # The accession gate fires for either label; the positive gate only for label 1.
    assert _by_name(results, "SaPS: train/test accession overlap").verdict == "fail"


def test_leakage_allows_cross_task_sharing(tmp_path):
    tables = _tables()
    tables[("base", "train")].append(
        {"UniprotEntry": "SHARE1", "task": "SaPS", "split": "train",
         "source_sheet": "NoPS", "label": "0"}
    )
    tables[("base", "train")].append(
        {"UniprotEntry": "SHARE1", "task": "PdPS", "split": "train",
         "source_sheet": "NoPS", "label": "0"}
    )
    _write_all(tmp_path, tables)
    results = leakage_layer(root=tmp_path)
    assert not [r for r in results if r.hard and r.verdict == "fail"]
    row = _by_name(results, "base: cross-task train accession sharing")
    assert row.verdict == "pass"
    assert row.hard is False
    assert row.measured["count"] >= 1


def test_leakage_allows_shared_test_sheets(tmp_path):
    tables = _tables()
    tables[("base", "test")].append(
        {"UniprotEntry": "TSHARE1", "task": "SaPS", "split": "test",
         "source_sheet": "PS-test", "label": "1"}
    )
    tables[("base", "test")].append(
        {"UniprotEntry": "TSHARE1", "task": "PdPS", "split": "test",
         "source_sheet": "PS-test", "label": "1"}
    )
    _write_all(tmp_path, tables)
    results = leakage_layer(root=tmp_path)
    assert not [r for r in results if r.hard and r.verdict == "fail"]
    row = _by_name(results, "base: cross-task test accession sharing")
    assert row.verdict == "pass"
    assert row.hard is False
    assert row.measured["count"] >= 1
    assert "PS-test" in row.measured["sheets"]


def test_leakage_sequence_check_skips_without_sequences(tmp_path):
    _write_all(tmp_path, _tables())
    row = _by_name(
        leakage_layer(root=tmp_path, sequences=None), "SaPS: train/test sequence overlap"
    )
    assert row.verdict == "skip"
    assert row.verdict != "pass"
    assert "uniprot_cache" in row.attribution


def test_leakage_sequence_identity_is_reported_without_blocking(tmp_path):
    tables = _tables()
    tables[("base", "train")].append(
        {"UniprotEntry": "SEQTR", "task": "SaPS", "split": "train",
         "source_sheet": "SaPS", "label": "1"}
    )
    tables[("base", "test")].append(
        {"UniprotEntry": "SEQTE", "task": "SaPS", "split": "test",
         "source_sheet": "SaPS-test", "label": "1"}
    )
    _write_all(tmp_path, tables)
    sequences = {"SEQTR": "MSHARED", "SEQTE": "MSHARED"}
    row = _by_name(
        leakage_layer(root=tmp_path, sequences=sequences), "SaPS: train/test sequence overlap"
    )
    assert row.verdict == "fail"
    assert row.hard is False
    assert row.measured["shared_sequences"] == 1
    assert row.measured["accession_pairs"] == 1
    assert row.measured["label_conflict_sequences"] == 0
    example = row.measured["examples"][0]
    assert example["train_accession"] == "SEQTR"
    assert example["test_accession"] == "SEQTE"
    assert example["train_label"] == [1]
    assert example["test_label"] == [1]
    assert "orthologue" in row.attribution
    assert "S2/S3" in row.attribution
    assert exit_code_for([row]) == 0
    assert exit_code_for([row], strict_leakage=True) == 1


def test_leakage_sequence_check_reports_the_label_conflict_subset(tmp_path):
    tables = _tables()
    tables[("base", "train")].append(
        {"UniprotEntry": "POS_TR", "task": "SaPS", "split": "train",
         "source_sheet": "SaPS", "label": "1"}
    )
    tables[("base", "test")].append(
        {"UniprotEntry": "NEG_TE", "task": "SaPS", "split": "test",
         "source_sheet": "NoPS-test", "label": "0"}
    )
    _write_all(tmp_path, tables)
    sequences = {"POS_TR": "MCONFLICT", "NEG_TE": "MCONFLICT"}
    row = _by_name(
        leakage_layer(root=tmp_path, sequences=sequences), "SaPS: train/test sequence overlap"
    )
    assert row.verdict == "fail"
    assert row.hard is False
    assert row.measured["shared_sequences"] == 1
    assert row.measured["accession_pairs"] == 1
    assert row.measured["label_conflict_sequences"] == 1
    conflict = row.measured["label_conflict_examples"][0]
    assert conflict["train_accession"] == "POS_TR"
    assert conflict["train_label"] == [1]
    assert conflict["test_accession"] == "NEG_TE"
    assert conflict["test_label"] == [0]
    assert "test AUC" in row.attribution


# ------------------------------------------------------ consistency layer


def test_consistency_flags_a_wrong_feature_definitions_version(tmp_path):
    _make_all(tmp_path)
    _make_artifact(tmp_path, "SaPS", feature_definitions_version="v2021")
    results = consistency_layer(models_dir=tmp_path)
    assert any(
        r.hard and r.verdict == "fail" and "feature_definitions_version" in r.name
        for r in results
    )


def test_consistency_flags_a_wrong_model_count(tmp_path):
    _make_all(tmp_path)
    _make_artifact(tmp_path, "SaPS", n_files=9)
    results = consistency_layer(models_dir=tmp_path)
    assert any(r.hard and r.verdict == "fail" and "model_count" in r.name for r in results)


def test_consistency_flags_a_wrong_seed(tmp_path):
    _make_all(tmp_path)
    _make_artifact(tmp_path, "SaPS", seed=7)
    results = consistency_layer(models_dir=tmp_path)
    assert any(r.hard and r.verdict == "fail" and r.name.endswith(": seed") for r in results)


def test_consistency_reports_a_version_drift_as_a_soft_note(tmp_path):
    _make_all(tmp_path)
    _make_artifact(tmp_path, "SaPS", xgboost_version="0.1.0")
    results = consistency_layer(models_dir=tmp_path)
    row = _by_name(results, "SaPS: xgboost_version")
    assert row.verdict == "fail"
    assert row.hard is False
    assert exit_code_for([row]) == 0
    assert exit_code_for([row], strict_paper=True) == 0


def test_consistency_flags_a_missing_manifest(tmp_path):
    _make_all(tmp_path)
    (tmp_path / "SaPS" / "manifest.json").unlink()
    results = consistency_layer(models_dir=tmp_path)
    failures = [r for r in results if r.hard and r.verdict == "fail"]
    assert failures
    assert any(
        "manifest" in r.attribution.lower() and "regenerate" in r.attribution.lower()
        for r in failures
    )


def test_consistency_passes_on_a_well_formed_artifact(tmp_path):
    _make_all(tmp_path)
    results = consistency_layer(models_dir=tmp_path)
    assert results
    assert all(r.verdict == "pass" for r in results), [r.name for r in results]


# ------------------------------------------------------- report and exit


def test_exit_code_is_zero_when_only_paper_rows_fail():
    results = [_result(verdict="fail", hard=False)]
    assert exit_code_for(results) == 0


def test_exit_code_is_one_on_a_hard_failure():
    results = [_result(layer="leakage", verdict="fail", hard=True)]
    assert exit_code_for(results) == 1


def test_strict_paper_promotes_paper_failures():
    results = [_result(verdict="fail", hard=False)]
    assert exit_code_for(results) == 0
    assert exit_code_for(results, strict_paper=True) == 1


def test_strict_paper_promotes_only_the_paper_layer():
    note = _result(
        layer="consistency",
        name="SaPS: xgboost_version",
        verdict="fail",
        hard=False,
    )
    assert exit_code_for([note], strict_paper=True) == 0
    paper = _result(layer="paper", name="hPdPS-10", verdict="fail", hard=False)
    assert exit_code_for([paper], strict_paper=True) == 1


def test_exit_code_ignores_a_sequence_identity_failure_without_strict_leakage():
    row = _result(
        layer="leakage",
        name="SaPS: train/test sequence overlap",
        verdict="fail",
        hard=False,
    )
    assert exit_code_for([row]) == 0
    assert exit_code_for([row], strict_leakage=True) == 1


def test_strict_flags_do_not_promote_a_consistency_note():
    note = _result(
        layer="consistency",
        name="SaPS: xgboost_version",
        verdict="fail",
        hard=False,
    )
    assert exit_code_for([note], strict_paper=True) == 0
    assert exit_code_for([note], strict_leakage=True) == 0


def test_relative_never_emits_an_absolute_path(tmp_path):
    outside = tmp_path / "outside_metrics.json"
    rendered = _relative(outside)
    assert not rendered.startswith("/")
    assert "outside_metrics.json" in rendered
    assert "outside repository" in rendered


def test_exit_code_is_zero_when_a_hard_check_skips():
    results = [_result(layer="leakage", verdict="skip", hard=True)]
    assert exit_code_for(results) == 0


def test_report_shape_is_stable():
    report = results_to_report([_result(measured=0.8, threshold=0.9)])
    assert set(report) == {"generated_from", "checks", "summary", "thresholds", "inputs"}
    check_keys = {"layer", "name", "verdict", "hard", "measured", "threshold", "attribution"}
    assert report["checks"]
    for check in report["checks"]:
        assert set(check) == check_keys
    assert set(report["summary"]) == {"pass", "fail", "skip", "hard_fail"}
    assert not _walk_keys(report) & {"generated_at", "timestamp"}


def test_format_report_table_names_every_failing_row():
    results = [
        _result(name="SaPS-8", verdict="pass", measured=0.861, threshold=0.862),
        _result(layer="leakage", name="hPdPS: leak", verdict="fail", hard=True,
                measured="acc123", threshold=0),
        _result(name="hPdPS-10", verdict="fail", hard=False, measured=0.8162,
                threshold=0.827, attribution="per-model; S5"),
    ]
    text = format_report_table(results)
    assert "hPdPS: leak" in text
    assert "acc123" in text
    assert "hPdPS-10" in text
    assert "0.8162" in text
    assert "SaPS-8" in text


# --------------------------------------------------------------- CLI driver

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_phasepred.py"


def _load_driver():
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location("validate_phasepred", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["validate_phasepred"] = module
    spec.loader.exec_module(module)
    return module


def test_cli_exits_two_without_a_traceback_on_a_missing_table(tmp_path, capsys):
    driver = _load_driver()
    empty = tmp_path / "empty"
    empty.mkdir()
    output = tmp_path / "out.json"
    code = driver.main(["--root", str(empty), "--output", str(output)])
    captured = capsys.readouterr()
    assert code == 2
    assert "Traceback" not in captured.err
    assert "chen2022_s2s3_base-features_train.tsv" in captured.err
    assert not output.exists()
