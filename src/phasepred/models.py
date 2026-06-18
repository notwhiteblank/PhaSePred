from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
