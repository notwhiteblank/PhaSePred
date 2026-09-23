"""Distilled XGBoost surrogate scorer (``mode="distilled"`` artifacts).

The surrogate is a GBDT regressor (locked S4 configuration, see
``provenance`` in the distilled artifact) that predicts the PhaSePred
web-archive catGRANULE ``single`` from the rich 181-dim sequence feature map
in :mod:`catgranule._distill`. The model is a pure-sequence computation with
no external tools; it ships inside the package (``weights/*.ubj``). Scoring
uses only ``xgboost`` (Booster) and never requires scikit-learn.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xgboost as xgb

from catgranule._distill import featurize
from catgranule.weights import CatGranuleWeights

_surrogate_cache: dict[Path, DistilledSurrogate] = {}

FEATURE_DIM = 181


class DistilledSurrogateError(RuntimeError):
    """Raised when a surrogate cannot be built or used."""


class DistilledSurrogate:
    """Lazy XGBoost surrogate bound to a ``mode="distilled"`` artifact."""

    def __init__(self, weights: CatGranuleWeights) -> None:
        if not weights.is_distilled:
            raise DistilledSurrogateError(
                "DistilledSurrogate requires a mode='distilled' artifact "
                f"(got route={weights.provenance.route!r}, mode={weights.mode!r})"
            )
        if weights.model_path is None:
            raise DistilledSurrogateError("Distilled artifact is missing its model file")
        self.weights = weights
        self._booster: xgb.Booster | None = None

    @property
    def model(self) -> xgb.Booster:
        if self._booster is None:
            booster = xgb.Booster()
            booster.load_model(str(self.weights.model_path))
            self._booster = booster
        return self._booster

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict an (N, 181) feature array -> (N,) float scores."""
        if X.shape[1] != FEATURE_DIM:
            raise DistilledSurrogateError(
                f"Expected {FEATURE_DIM} features, got {X.shape[1]}"
            )
        return self.model.predict(xgb.DMatrix(X))

    def predict_single(self, sequence: str) -> float:
        vec = featurize(sequence, self.weights).reshape(1, FEATURE_DIM)
        return float(self.model.predict(xgb.DMatrix(vec))[0])

    @property
    def model_path(self) -> Path:
        assert self.weights.model_path is not None
        return self.weights.model_path


def load_surrogate(weights: CatGranuleWeights) -> DistilledSurrogate:
    """Return a process-wide cached surrogate for the given distilled artifact."""
    if weights.model_path is None:
        raise DistilledSurrogateError("Distilled artifact is missing its model file")
    cached = _surrogate_cache.get(weights.model_path)
    if cached is None:
        cached = DistilledSurrogate(weights)
        _surrogate_cache[weights.model_path] = cached
    return cached


__all__ = [
    "DistilledSurrogate",
    "DistilledSurrogateError",
    "load_surrogate",
    "FEATURE_DIM",
]