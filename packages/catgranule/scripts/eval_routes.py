"""S4 gate evaluation for the three candidate routes (CATGRANULE_SPEC 4.2-4.3).

Routes (all evaluated on the web-archive human holdout, seed 42):

  --route paper     : formula reproduction with the published coefficients and a
                      Z-normalization recomputed from the local yeast proteome
                      (SPEC 4.1 main route).
  --route linear7   : the "scaled weights" control proposed in SPEC 4.1 - OLS on the
                      [mean of each of the six heptapeptide scale means, log L]
                      feature vector, trained on the web eval-pool. Reported as a
                      comparison, not a product.
  --route distilled : the distilled model (SPEC 4.3, D11 pre-authorized): XGBoost
                      regression on the rich feature map, trained on the eval-pool
                      only. Evaluated exactly once on the holdout.

Gate (PLAN sec.5): Spearman rho >= 0.99 AND median |delta| <= 0.01 (single).
Residue-level consistency is reported separately (see CATGRANULE_VALIDATION.md).

Run:
  uv run python packages/catgranule/scripts/eval_routes.py --route paper
  uv run python packages/catgranule/scripts/eval_routes.py --route linear7
  uv run python packages/catgranule/scripts/eval_routes.py --route distilled

Outputs are written to data/raw/external/phasepred_web/ (gitignored).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.stats as st
import xgboost as xgb
from catgranule.weights import load_weights
from sklearn.linear_model import LinearRegression

_scripts_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(_scripts_dir))
from _s4_numeric import (  # noqa: E402
    _STANDARD_AA,
    encode_bytes,
    featurize_distill,
    make_scale_matrix,
    paper_coefficients,
    raw_score,
    window_mean,
)

WEB_DIR = Path("data/raw/external/phasepred_web")
CACHE = Path("data/interim/uniprot_cache.jsonl")
NORM_PAPER = WEB_DIR / "paper_route_yeast_norm.json"
PAPER_HOLDOUT = WEB_DIR / "paper_route_holdout_results.json"
DISTILL_HOLDOUT = WEB_DIR / "distilled_holdout_results.json"


def _read_compact(species: str, want: set[str] | None = None) -> list[dict]:
    rows = []
    with (WEB_DIR / f"{species}_compact.jsonl").open() as fh:
        for line in fh:
            r = json.loads(line)
            if want is not None and r["accession"] not in want:
                continue
            rows.append(r)
    return rows


def yeast_norm(scale_matrix: np.ndarray, coefficients: np.ndarray) -> dict:
    """SPEC 4.1: recompute Z constants over the local yeast proteome (cache)."""
    raws: list[float] = []
    skipped = 0
    with CACHE.open() as fh:
        for line in fh:
            rec = json.loads(line)
            if "Saccharomyces cerevisiae" not in rec["header"]:
                continue
            seq = rec["sequence"]
            if len(seq) == 0 or not set(seq) <= _STANDARD_AA:
                skipped += 1
                continue
            raws.append(raw_score(seq, scale_matrix, coefficients))
    raws = np.array(raws)
    norm = {
        "mean": float(raws.mean()),
        "std": float(raws.std()),
        "n_proteins": int(len(raws)),
        "n_skipped": int(skipped),
        "source": "data/interim/uniprot_cache.jsonl (all OS=Saccharomyces cerevisiae)",
        "method": "Bolognesi 2016 Equation 1-2 with published coefficients; "
        "mean/std of raw Equation-2 values over the yeast proteome",
    }
    NORM_PAPER.write_text(json.dumps(norm, indent=2))
    print(json.dumps(norm, indent=2), flush=True)
    return norm


def paper_route() -> None:
    weights = load_weights()
    scale_matrix = make_scale_matrix(weights)
    coefficients = paper_coefficients(weights)
    norm = yeast_norm(scale_matrix, coefficients)
    mu, sd = norm["mean"], norm["std"]

    split = json.loads((WEB_DIR / "holdout_human.json").read_text())
    hold = set(split["holdout"])
    rows = []
    for r in _read_compact("human", hold):
        seq = r["sequence"]
        if not set(seq) <= _STANDARD_AA:
            continue
        raw = raw_score(seq, scale_matrix, coefficients)
        rows.append((r["accession"], (raw - mu) / sd, r["cg_single"]))
    pred = np.array([x[1] for x in rows])
    web = np.array([x[2] for x in rows])
    rho = float(st.spearmanr(pred, web).correlation)
    pearson = float(np.corrcoef(pred, web)[0, 1])
    delta = pred - web
    ad = np.abs(delta)
    result = {
        "route": "paper",
        "n_holdout": int(len(rows)),
        "spearman_rho": round(rho, 6),
        "pearson_r": round(pearson, 6),
        "median_abs_delta": round(float(np.median(ad)), 6),
        "mean_abs_delta": round(float(ad.mean()), 6),
        "q90_abs_delta": round(float(np.quantile(ad, 0.90)), 6),
        "gate_rho_pass": bool(rho >= 0.99),
        "gate_delta_pass": bool(float(np.median(ad)) <= 0.01),
        "gate_pass": bool(rho >= 0.99 and float(np.median(ad)) <= 0.01),
        "normalization": norm,
        "coefficients": "published (aR=(0.48,7.24) aD=(0.26,11.54) aP=(1.98,1.42) aL=0.25)",
        "seed": 42,
    }
    PAPER_HOLDOUT.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


def linear7_control() -> None:
    """SPEC 4.1 control: scaled-weights linear refit on the 7-vector."""
    weights = load_weights()
    scale_matrix = make_scale_matrix(weights)
    split = json.loads((WEB_DIR / "holdout_human.json").read_text())
    evl, hold = set(split["eval_pool"]), set(split["holdout"])

    def vec7(seq: str) -> list[float]:
        vals = scale_matrix[encode_bytes(seq)]
        means = window_mean(vals, 3).mean(axis=0)
        return list(means) + [math.log(len(seq))]

    rows = []
    for want in (evl, hold):
        for r in _read_compact("human", want):
            if not set(r["sequence"]) <= _STANDARD_AA:
                continue
            rows.append((tuple(vec7(r["sequence"])), r["cg_single"], r["accession"]))
    arr = np.array([x[0] for x in rows])
    y = np.array([x[1] for x in rows])
    accs = np.array([x[2] for x in rows])
    is_eval = np.array([a in evl for a in accs])
    model = LinearRegression().fit(arr[is_eval], y[is_eval])
    pred = model.predict(arr[~is_eval])
    web = y[~is_eval]
    delta = pred - web
    ad = np.abs(delta)
    result = {
        "route": "linear7-control",
        "n_holdout": int(len(web)),
        "spearman_rho": round(float(st.spearmanr(pred, web).correlation), 6),
        "pearson_r": round(float(np.corrcoef(pred, web)[0, 1]), 6),
        "median_abs_delta": round(float(np.median(ad)), 6),
        "q90_abs_delta": round(float(np.quantile(ad, 0.90)), 6),
        "fitted_coefficients": np.round(model.coef_, 4).tolist(),
        "intercept": round(float(model.intercept_), 4),
        "note": "control only; trained on eval-pool (web single target), not a product",
    }
    out = WEB_DIR / "linear7_holdout_results.json"
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


def distilled_route() -> None:
    """SPEC 4.3 distilled model: XGBoost on the rich feature map (locked config)."""
    weights = load_weights()
    scale_matrix = make_scale_matrix(weights)
    coefficients = paper_coefficients(weights)
    split = json.loads((WEB_DIR / "holdout_human.json").read_text())
    evl, hold = set(split["eval_pool"]), set(split["holdout"])
    want = evl | hold
    accs: list[str] = []
    X: list[np.ndarray] = []
    y: list[float] = []
    t0 = time.time()
    for r in _read_compact("human", want):
        seq = r["sequence"]
        if len(seq) < 60 or not set(seq) <= _STANDARD_AA:
            continue
        accs.append(r["accession"])
        X.append(featurize_distill(seq, scale_matrix, coefficients))
        y.append(float(r["cg_single"]))
    X = np.nan_to_num(np.array(X), nan=0.0, posinf=0.0, neginf=0.0)
    y = np.array(y)
    print(f"featurized {len(accs)} proteins in {time.time()-t0:.1f}s", flush=True)
    accs = np.array(accs)
    is_eval = np.array([a in evl for a in accs])
    Xe, ye = X[is_eval], y[is_eval]
    Xh, yh = X[~is_eval], y[~is_eval]
    print("eval", Xe.shape, "holdout", Xh.shape, flush=True)

    model = xgb.XGBRegressor(
        n_estimators=600,
        max_depth=7,
        learning_rate=0.04,
        subsample=0.85,
        colsample_bytree=0.8,
        min_child_weight=3,
        reg_lambda=5.0,
        n_jobs=64,
        random_state=0,
    )
    model.fit(Xe, ye)
    pred = model.predict(Xh)
    rho = float(st.spearmanr(pred, yh).correlation)
    pearson = float(np.corrcoef(pred, yh)[0, 1])
    delta = yh - pred
    ad = np.abs(delta)
    result = {
        "route": "distilled",
        "model": "XGBRegressor(n_estimators=600, max_depth=7, lr=0.04, "
        "subsample=0.85, colsample_bytree=0.8, random_state=0)",
        "features": "rich map: window summaries radii 3/7/13/25, g-profile stats "
        "half-windows 3/12/25, AA composition, L/logL",
        "trained_on": "web-archive human eval-pool only",
        "train_n": int(is_eval.sum()),
        "n_holdout": int(len(yh)),
        "spearman_rho": round(rho, 6),
        "pearson_r": round(pearson, 6),
        "r2": round(float(1 - delta.var() / yh.var()), 6),
        "median_abs_delta": round(float(np.median(ad)), 6),
        "mean_abs_delta": round(float(ad.mean()), 6),
        "q90_abs_delta": round(float(np.quantile(ad, 0.90)), 6),
        "gate_rho_pass": bool(rho >= 0.99),
        "gate_delta_pass": bool(float(np.median(ad)) <= 0.01),
        "gate_pass": bool(rho >= 0.99 and float(np.median(ad)) <= 0.01),
        "seed_fixed": True,
    }
    DISTILL_HOLDOUT.write_text(json.dumps(result, indent=2))
    np.savez(
        WEB_DIR / "distilled_holdout_predictions.npz",
        acc=accs[~is_eval],
        web=yh,
        pred=pred,
    )
    print(json.dumps(result, indent=2), flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--route", choices=["paper", "linear7", "distilled"], default="paper")
    args = ap.parse_args()
    if args.route == "paper":
        paper_route()
    elif args.route == "linear7":
        linear7_control()
    else:
        distilled_route()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
