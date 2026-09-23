"""PhaSePred predictor: compute features and score proteins from FASTA."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from phasepred.data import models_root

# Default artifact (single product): Product A — paper split recomputed.
DEFAULT_MODELS_DIR = models_root()


# ---------------------------------------------------------------------------
# Feature computation
# ---------------------------------------------------------------------------


def compute_all_features(
    records: list[dict[str, str]],
    *,
    include_human: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """Compute all PhaSePred features for a list of FASTA records.

    Each record is {"accession": str, "sequence": str}.
    Returns a DataFrame with feature columns ready for prediction.
    """
    from catgranule import CatGranuleInputError, score_sequence

    from phasepred import tools
    from phasepred.features import compute_native_features as _native

    t_start = time.time()

    # --- Native features (pure Python, fast) ---
    native_rows = []
    for r in records:
        seq = str(r["sequence"]).upper()
        try:
            feats = _native(seq)
            native_rows.append(
                {
                    "UniprotEntry": r["accession"],
                    "length": feats["length"],
                    "Hydropathy": feats["Hydropathy"],
                    "FCR": feats["FCR"],
                }
            )
        except Exception:
            native_rows.append(
                {
                    "UniprotEntry": r["accession"],
                    "length": float("nan"),
                    "Hydropathy": float("nan"),
                    "FCR": float("nan"),
                }
            )

    feature_df = pd.DataFrame(native_rows)
    if verbose:
        print(f"  Native features: {len(feature_df)} in {time.time() - t_start:.1f}s")

    # --- catGRANULE (via the catgranule package, legacy-unaudited weights) ---
    t0 = time.time()
    cat_scores: dict[str, float] = {}
    for r in records:
        try:
            cat_scores[r["accession"]] = float(score_sequence(str(r["sequence"]))["single"])
        except CatGranuleInputError:
            cat_scores[r["accession"]] = float("nan")
    feature_df["catGRANULE"] = feature_df["UniprotEntry"].map(cat_scores)
    if verbose:
        n = feature_df["catGRANULE"].notna().sum()
        print(f"  catGRANULE: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- ESpritz IDR ---
    t0 = time.time()
    espritz_scores = tools.run_espritz(records)
    feature_df["IDR"] = feature_df["UniprotEntry"].map(espritz_scores)
    if verbose:
        n = feature_df["IDR"].notna().sum()
        print(f"  ESpritz IDR: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- SEG LCR ---
    t0 = time.time()
    seg_scores = tools.run_seg(records)
    feature_df["LCR"] = feature_df["UniprotEntry"].map(seg_scores)
    if verbose:
        n = feature_df["LCR"].notna().sum()
        print(f"  SEG LCR: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- PScore (parallel) ---
    t0 = time.time()
    pscore_scores = tools.run_pscore(records)
    feature_df["PScore"] = feature_df["UniprotEntry"].map(pscore_scores)
    if verbose:
        n = feature_df["PScore"].notna().sum()
        print(f"  PScore: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- PLAAC (parallel) ---
    t0 = time.time()
    plaac_scores = tools.run_plaac(records)
    feature_df["PLAAC"] = feature_df["UniprotEntry"].map(plaac_scores)
    if verbose:
        n = feature_df["PLAAC"].notna().sum()
        print(f"  PLAAC: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- DeepCoil (on-demand subprocess into isolated conda env; binary 0/1) ---
    t0 = time.time()
    dc_scores = tools.run_deepcoil(records)
    feature_df["DeepCoil"] = feature_df["UniprotEntry"].map(dc_scores)
    if verbose:
        n = feature_df["DeepCoil"].notna().sum()
        print(f"  DeepCoil: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- Human-specific features ---
    if include_human:
        t0 = time.time()
        phos_scores = tools.compute_phos_freq(records)
        feature_df["Phos freq"] = feature_df["UniprotEntry"].map(phos_scores)
        dp_scores = tools.lookup_deepphase([r["accession"] for r in records])
        feature_df["DeepPhase"] = feature_df["UniprotEntry"].map(dp_scores)
        if verbose:
            n_p = feature_df["Phos freq"].notna().sum()
            n_d = feature_df["DeepPhase"].notna().sum()
            print(f"  Phos freq + DeepPhase: {n_p}, {n_d} in {time.time() - t0:.1f}s")

    if verbose:
        total_t = time.time() - t_start
        print(f"  Total feature computation: {total_t:.0f}s")

    return feature_df


# ---------------------------------------------------------------------------
# Model loading and prediction
# ---------------------------------------------------------------------------


def load_ensemble(model_dir: Path) -> list:
    """Load a 10-model XGBoost ensemble from a directory of joblib files."""
    import joblib

    model_files = sorted(model_dir.glob("8f_model_*.joblib"))
    if not model_files:
        raise FileNotFoundError(f"No 8f_model_*.joblib files found in {model_dir}")
    return [joblib.load(str(f)) for f in model_files]


def load_single_model(model_path: Path):
    """Load a single model artifact (joblib file)."""
    import joblib

    if not model_path.exists():
        raise FileNotFoundError(f"Model not found: {model_path}")
    return joblib.load(str(model_path))


def predict_scores(
    models: list,
    features: pd.DataFrame,
    feature_columns: list[str],
) -> pd.DataFrame:
    """Predict scores using a loaded ensemble (or single model).

    If the model was trained on a subset of columns (e.g. paper baseline
    with only 3 features), use only those columns from the feature table.

    Returns DataFrame with columns: UniprotEntry, score
    """
    # Check if the model artifact specifies its own feature columns
    model_columns = None
    if isinstance(models, list) and models:
        first = models[0]
        if hasattr(first, "feature_columns"):
            model_columns = list(first.feature_columns)
    elif hasattr(models, "feature_columns"):
        model_columns = list(models.feature_columns)

    effective_columns = model_columns or feature_columns
    avail_cols = [c for c in effective_columns if c in features.columns]
    missing = set(effective_columns) - set(avail_cols)
    if missing:
        for col in missing:
            features[col] = float("nan")

    X = features[avail_cols]

    # models may be a single object or a list
    model_list = models if isinstance(models, list) else [models]

    if len(model_list) > 1:
        scores = np.zeros(len(X))
        for model in model_list:
            scores += _model_predict_proba(model, X)[:, 1]
        scores /= len(model_list)
    else:
        scores = _model_predict_proba(model_list[0], X)[:, 1]

    result = features[["UniprotEntry"]].copy()
    result["score"] = scores
    return result


def _model_predict_proba(model, X: pd.DataFrame) -> np.ndarray:
    """Call predict_proba, handling both sklearn Pipeline and ModelArtifact."""
    if hasattr(model, "estimator"):
        # ModelArtifact wrapper
        return model.estimator.predict_proba(X)
    return model.predict_proba(X)