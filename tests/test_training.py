"""Protocol-level unit tests for the one-command retrain (E2, Task 1).

These cover the pieces of ``phasepred.training`` that do not require a full
training run: negative-set sampling, feature-column scope, the committed
tables' shape, the positive/negative split, leakage detection, ensemble
scoring, and the paper-default XGBoost constructor. Training the real
ensembles is gated by G1 rather than by pytest.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from xgboost import XGBClassifier

from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS
from phasepred.training import (
    N_MODELS,
    NEG_RATIO,
    check_leakage,
    ensemble_score,
    feature_columns,
    make_xgb,
    neg_sets,
    split_pos_neg,
    task_frame,
)


def _pool(size: int = 50) -> pd.DataFrame:
    return pd.DataFrame(
        {"UniprotEntry": [f"P{i:05d}" for i in range(size)], "value": range(size)}
    )


def test_neg_sets_are_deterministic():
    pool = _pool()
    first = neg_sets(pool, 10, 42)
    second = neg_sets(pool, 10, 42)
    assert len(first) == N_MODELS
    for left, right in zip(first, second, strict=True):
        assert left["UniprotEntry"].tolist() == right["UniprotEntry"].tolist()


def test_neg_sets_differ_per_index():
    pool = _pool()
    sets = neg_sets(pool, 10, 42)
    expected_size = min(len(pool), 10 * NEG_RATIO)
    assert all(len(neg) == expected_size for neg in sets)
    signatures = {tuple(neg["UniprotEntry"].tolist()) for neg in sets}
    assert len(signatures) == N_MODELS


def test_neg_sets_honour_exclude():
    pool = _pool()
    excluded_index = pool.index[:20]
    sets = neg_sets(pool, 10, 42, exclude=excluded_index)
    excluded = set(pool.loc[excluded_index, "UniprotEntry"])
    for neg in sets:
        assert not set(neg["UniprotEntry"]) & excluded


def test_feature_columns_follow_the_scope():
    assert len(BASE_FEATURE_COLUMNS) == 8
    assert len(HUMAN_FEATURE_COLUMNS) == 10
    assert feature_columns("SaPS") == list(BASE_FEATURE_COLUMNS)
    assert feature_columns("PdPS") == list(BASE_FEATURE_COLUMNS)
    assert feature_columns("hSaPS") == list(HUMAN_FEATURE_COLUMNS)
    assert feature_columns("hPdPS") == list(HUMAN_FEATURE_COLUMNS)


def test_task_frame_reads_the_committed_tables():
    # Row counts are metrics.json's train_rows/test_rows. Positive counts are
    # measured from the tables and are NOT derivable from metrics.json's
    # n_pos_train: test-side positives come from merging the *-test sheet with
    # PS-test / hPS-test, so they differ from the train-side counts.
    # Phase E3's leakage check reuses these eight (task, split) pairs.
    expected = {
        ("SaPS", "train"): (48286, 128),
        ("SaPS", "test"): (12188, 126),
        ("PdPS", "train"): (48372, 214),
        ("PdPS", "test"): (12228, 166),
        ("hSaPS", "train"): (8860, 59),
        ("hSaPS", "test"): (2257, 57),
        ("hPdPS", "train"): (8897, 96),
        ("hPdPS", "test"): (2283, 83),
    }
    for (task, split), (rows, positives) in expected.items():
        frame = task_frame(task, split)
        assert len(frame) == rows, (task, split)
        assert int((frame["label"] == 1).sum()) == positives, (task, split)


def test_split_pos_neg_casts_to_float():
    frame = task_frame("SaPS", "train")
    columns = feature_columns("SaPS")
    pos, neg = split_pos_neg(frame, columns)
    assert len(pos) == int((frame["label"] == 1).sum())
    assert len(neg) == int((frame["label"] == 0).sum())
    for out in (pos, neg):
        assert list(out.columns) == columns
        for column in columns:
            assert out[column].dtype == np.float64, column


def test_check_leakage_counts_overlap():
    train = pd.DataFrame({"UniprotEntry": ["A", "B", "C"]})
    test = pd.DataFrame({"UniprotEntry": ["C", "D"]})
    assert check_leakage(train, test) == 1
    assert check_leakage(train, pd.DataFrame({"UniprotEntry": ["X", "Y"]})) == 0


def test_ensemble_score_is_the_mean_of_proba():
    rng = np.random.default_rng(0)
    features = pd.DataFrame(rng.normal(size=(30, 3)), columns=["a", "b", "c"])
    labels = np.array([0, 1] * 15)
    models = [make_xgb(42), make_xgb(42)]
    for model in models:
        model.fit(features, labels)
    probs = np.column_stack([model.predict_proba(features)[:, 1] for model in models])
    np.testing.assert_array_equal(ensemble_score(models, features), probs.mean(axis=1))


def test_make_xgb_uses_paper_defaults():
    model = make_xgb(42)
    assert model.random_state == 42
    assert model.n_jobs == 4
    assert model.eval_metric == "logloss"
    default = XGBClassifier()
    for name in ("scale_pos_weight", "max_depth", "learning_rate", "n_estimators", "objective"):
        assert model.get_params()[name] == default.get_params()[name], name
