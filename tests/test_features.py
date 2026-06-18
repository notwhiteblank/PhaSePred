from __future__ import annotations

import math

import pandas as pd
import pytest

from phasepred.features import (
    BASE_FEATURE_COLUMNS,
    HUMAN_FEATURE_COLUMNS,
    FeatureSchemaError,
    build_feature_matrix,
    compute_native_features,
)


def test_compute_native_features_returns_length_hydropathy_and_fcr() -> None:
    features = compute_native_features("ACDE")

    assert features["length"] == 4
    assert features["FCR"] == 0.5
    assert math.isclose(features["Hydropathy"], 0.425, rel_tol=1e-9)


def test_build_feature_matrix_joins_required_columns() -> None:
    records = pd.DataFrame(
        {
            "UniprotEntry": ["P1", "P2"],
            "sequence": ["ACDE", "KKKK"],
        }
    )
    imported = pd.DataFrame(
        {
            "UniprotEntry": ["P1", "P2"],
            "IDR": [0.1, 0.2],
            "LCR": [0.3, 0.4],
            "PScore": [1.0, 2.0],
            "PLAAC": [0.5, 0.6],
            "catGRANULE": [0.7, 0.8],
            "DeepCoil": [0.0, 1.0],
            "Phos freq": [0.9, 0.1],
            "DeepPhase": [0.2, 0.3],
        }
    )

    frame = build_feature_matrix(records, imported, feature_columns=HUMAN_FEATURE_COLUMNS)

    assert list(frame["UniprotEntry"]) == ["P1", "P2"]
    assert set(BASE_FEATURE_COLUMNS).issubset(frame.columns)
    assert set(HUMAN_FEATURE_COLUMNS).issubset(frame.columns)
    assert frame.loc[0, "length"] == 4


def test_build_feature_matrix_raises_on_missing_columns() -> None:
    records = pd.DataFrame(
        {
            "UniprotEntry": ["P1"],
            "sequence": ["ACDE"],
        }
    )
    imported = pd.DataFrame({"UniprotEntry": ["P1"], "IDR": [0.1]})

    with pytest.raises(FeatureSchemaError):
        build_feature_matrix(records, imported, feature_columns=BASE_FEATURE_COLUMNS)
