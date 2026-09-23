from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from phasepred.models import (
    FEATURE_DEFINITIONS_VERSION,
    ModelManifestError,
    load_model_artifact,
    load_model_manifest,
    model_dir_is_stale,
    predict_with_artifact,
    train_model_artifact,
)

pytestmark = pytest.mark.toolfree


def test_train_and_predict_round_trip(tmp_path: Path) -> None:
    features = pd.DataFrame(
        {
            "UniprotEntry": ["P1", "P2", "P3", "P4"],
            "Hydropathy": [1.0, 1.2, -1.0, -1.2],
            "FCR": [0.1, 0.2, 0.8, 0.9],
            "IDR": [0.0, 0.1, 0.8, 0.9],
            "LCR": [0.2, 0.2, 0.7, 0.8],
            "PScore": [0.1, 0.2, 0.9, 0.8],
            "PLAAC": [0.0, 0.1, 0.8, 0.9],
            "catGRANULE": [0.1, 0.1, 0.7, 0.8],
            "DeepCoil": [0.0, 0.1, 0.8, 0.9],
            "label": [0, 0, 1, 1],
        }
    )

    artifact = train_model_artifact(features, label_column="label", task="SaPS")
    model_path = tmp_path / "model.joblib"
    artifact.dump(model_path)
    loaded = load_model_artifact(model_path)

    predictions = predict_with_artifact(loaded, features.drop(columns=["label"]))

    assert list(predictions["UniprotEntry"]) == ["P1", "P2", "P3", "P4"]
    assert "score" in predictions.columns
    assert predictions["score"].between(0, 1).all()


def test_model_dir_without_manifest_is_stale(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()

    stale, reason = model_dir_is_stale(mode_dir)

    assert stale is True
    assert "manifest.json" in reason


def test_model_dir_fresh_manifest_not_stale(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()
    mode_dir.joinpath("manifest.json").write_text(
        json.dumps({"feature_definitions_version": FEATURE_DEFINITIONS_VERSION}),
        encoding="utf-8",
    )

    stale, reason = model_dir_is_stale(mode_dir)

    assert stale is False
    assert reason == ""


def test_model_dir_mismatched_version_is_stale(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()
    mode_dir.joinpath("manifest.json").write_text(
        json.dumps({"feature_definitions_version": "v2019"}),
        encoding="utf-8",
    )

    stale, reason = model_dir_is_stale(mode_dir)

    assert stale is True
    assert "v2019" in reason
    assert FEATURE_DEFINITIONS_VERSION in reason


def test_model_manifest_round_trip(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()
    manifest_path = mode_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps({"feature_definitions_version": FEATURE_DEFINITIONS_VERSION, "n_models": 10}),
        encoding="utf-8",
    )

    manifest = load_model_manifest(mode_dir)

    assert manifest is not None
    assert manifest.feature_definitions_version == FEATURE_DEFINITIONS_VERSION
    assert manifest.raw["n_models"] == 10


def test_model_manifest_missing_version_raises(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()
    mode_dir.joinpath("manifest.json").write_text(
        json.dumps({"n_models": 10}),
        encoding="utf-8",
    )

    with pytest.raises(ModelManifestError):
        load_model_manifest(mode_dir)


def test_model_manifest_invalid_json_raises(tmp_path: Path) -> None:
    mode_dir = tmp_path / "SaPS"
    mode_dir.mkdir()
    mode_dir.joinpath("manifest.json").write_text("not json{", encoding="utf-8")

    with pytest.raises(ModelManifestError):
        load_model_manifest(mode_dir)
