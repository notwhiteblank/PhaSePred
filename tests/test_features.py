from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from phasepred.features import (
    BASE_FEATURE_COLUMNS,
    HUMAN_FEATURE_COLUMNS,
    FeatureSchemaError,
    build_feature_matrix,
    compute_fcr,
    compute_native_features,
)

pytestmark = pytest.mark.toolfree


def test_compute_native_features_returns_length_hydropathy_and_fcr() -> None:
    features = compute_native_features("ACDE")

    assert features["length"] == 4
    assert features["FCR"] == 0.5
    assert math.isclose(features["Hydropathy"], 0.425, rel_tol=1e-9)


def test_fcr_is_defined_for_selenocysteine() -> None:
    assert compute_fcr("KUDER") == 0.8
    q8 = Path(__file__).parent / "fixtures" / "sequences" / "Q8WWX9.fasta"
    sequence = "".join(q8.read_text().splitlines()[1:])
    assert compute_fcr(sequence) == pytest.approx(0.23448275862068965)
    with pytest.raises(FeatureSchemaError, match="Unsupported residues"):
        compute_native_features("KUDER")


def test_predictor_keeps_fcr_when_hydropathy_is_unsupported(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from phasepred import tools
    from phasepred.predictor import compute_all_features

    for name in ("run_espritz", "run_seg", "run_pscore", "run_plaac", "run_deepcoil"):
        monkeypatch.setattr(tools, name, lambda records: {})
    row = compute_all_features(
        [{"accession": "SEL", "sequence": "KUDER"}], verbose=False
    ).iloc[0]
    assert row["length"] == 5
    assert row["FCR"] == 0.8
    assert pd.isna(row["Hydropathy"])
    assert pd.isna(row["catGRANULE"])


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
