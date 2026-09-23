from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import joblib
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from phasepred.features import (
    BASE_FEATURE_COLUMNS,
    HUMAN_FEATURE_COLUMNS,
    validate_feature_columns,
)

FEATURES_BY_TASK = {
    "SaPS": BASE_FEATURE_COLUMNS,
    "PdPS": BASE_FEATURE_COLUMNS,
    "hSaPS": HUMAN_FEATURE_COLUMNS,
    "hPdPS": HUMAN_FEATURE_COLUMNS,
}

# The v2022 feature-definitions version that the code implements
# (docs/FEATURES.md). Model artifacts tag their training definitions with
# this value; mismatches trigger the stale warning in `predict`.
FEATURE_DEFINITIONS_VERSION = "v2022"

MANIFEST_FILENAME = "manifest.json"


class ModelManifestError(RuntimeError):
    """Raised when a model directory manifest is unreadable."""


@dataclass(frozen=True)
class ModelManifest:
    feature_definitions_version: str
    raw: dict[str, Any]


def load_model_manifest(model_dir: Path) -> ModelManifest | None:
    """Read ``manifest.json`` from a model directory, if present.

    Returns ``None`` for pre-0.9.0 artifacts that carry no manifest. A
    malformed manifest raises :class:`ModelManifestError`.
    """
    manifest_path = model_dir / MANIFEST_FILENAME
    if not manifest_path.is_file():
        return None
    try:
        import json

        raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ModelManifestError(f"Unreadable model manifest at {manifest_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ModelManifestError(f"Model manifest at {manifest_path} is not a JSON object")
    version = raw.get("feature_definitions_version")
    if not isinstance(version, str) or not version:
        raise ModelManifestError(
            f"Model manifest at {manifest_path} is missing string "
            "feature_definitions_version"
        )
    return ModelManifest(feature_definitions_version=version, raw=raw)


def model_dir_is_stale(model_dir: Path) -> tuple[bool, str]:
    """Check whether a model directory binds an out-of-date feature set.

    Returns ``(is_stale, reason)``. Legacy artifacts without a manifest, and
    artifacts whose ``feature_definitions_version`` differs from the current
    ``FEATURE_DEFINITIONS_VERSION``, are stale.
    """
    manifest = load_model_manifest(model_dir)
    if manifest is None:
        return True, (
            f"model directory {model_dir} has no {MANIFEST_FILENAME} "
            "(pre-0.9.0 artifact trained before feature_definitions_version "
            f"tagging; current version {FEATURE_DEFINITIONS_VERSION})"
        )
    if manifest.feature_definitions_version != FEATURE_DEFINITIONS_VERSION:
        return True, (
            f"model directory {model_dir} binds feature_definitions_version "
            f"{manifest.feature_definitions_version!r} but the code implements "
            f"{FEATURE_DEFINITIONS_VERSION!r}"
        )
    return False, ""


@dataclass
class ModelArtifact:
    task: str
    feature_columns: list[str]
    estimator: Pipeline

    def dump(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path)


def train_model_artifact(
    frame: pd.DataFrame,
    *,
    label_column: str,
    task: str,
    feature_columns: list[str] | None = None,
) -> ModelArtifact:
    if label_column not in frame.columns:
        raise ValueError(f"Label column not found: {label_column}")
    columns = feature_columns or list(FEATURES_BY_TASK[task])
    validate_feature_columns(frame, columns)
    estimator = Pipeline(
        steps=[
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("classifier", LogisticRegression(max_iter=2000, random_state=42)),
        ]
    )
    estimator.fit(frame[columns], frame[label_column].astype(int))
    return ModelArtifact(task=task, feature_columns=columns, estimator=estimator)


def load_model_artifact(path: str | Path) -> ModelArtifact:
    artifact = joblib.load(path)
    if not isinstance(artifact, ModelArtifact):
        raise TypeError(f"Unexpected model artifact type: {type(artifact)!r}")
    return artifact


def predict_with_artifact(artifact: ModelArtifact, frame: pd.DataFrame) -> pd.DataFrame:
    validate_feature_columns(frame, artifact.feature_columns)
    scores = artifact.estimator.predict_proba(frame[artifact.feature_columns])[:, 1]
    output = pd.DataFrame({"score": scores})
    if "UniprotEntry" in frame.columns:
        output.insert(0, "UniprotEntry", frame["UniprotEntry"].to_numpy())
    return output