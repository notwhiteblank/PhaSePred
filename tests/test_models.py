from __future__ import annotations

from pathlib import Path

import pandas as pd

from phasepred.models import load_model_artifact, predict_with_artifact, train_model_artifact


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
