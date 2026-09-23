"""S4b: reproduce the distilled catGRANULE surrogate from the locked S4 config.

Locked configuration (docs/migration/work-log/S4.md breakpoint guide, D21):
  XGBRegressor(n_estimators=600, max_depth=7, learning_rate=0.04,
               subsample=0.85, colsample_bytree=0.8, min_child_weight=3,
               reg_lambda=5.0, random_state=0)
  Feature map: 181-dim (_s4_numeric.featurize_distill)
  Trained on web-archive human eval-pool only (n=16164); gate on holdout (n=4042).

Reproducibility gate (D21): retrain must give rho=0.998961 +/- 0.0005 and
median|delta|=0.017653 +/- 0.0005 on the holdout. If it does not match,
STOP (no artifacts written).

Writes the shipped package artifacts:
  src/catgranule/weights/distilled_xgb_model.ubj
  src/catgranule/weights/catgranule_distilled_weights.json
  src/catgranule/weights/catgranule_paper_weights.json
plus reports under data/raw/external/phasepred_web/ (gitignored).

Dependencies: the web-archive data area (S4 bootstrap) must be present
(data/raw/external/phasepred_web/human_compact.jsonl + holdout_human.json).

Run (from the repo root):
  uv run python packages/catgranule/scripts/reproduce_distilled.py --train-save
  uv run python packages/catgranule/scripts/reproduce_distilled.py --paper-artifact
  uv run python packages/catgranule/scripts/reproduce_distilled.py --cross-species
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.stats as st
import xgboost as xgb
from _s4_numeric import _STANDARD_AA, featurize_distill

REPO = Path(__file__).resolve().parents[3]
WEB_DIR = REPO / "data/raw/external/phasepred_web"
WEIGHTS_DIR = (
    REPO / "packages/catgranule/src/catgranule/weights"
)

LOCKED = dict(
    n_estimators=600,
    max_depth=7,
    learning_rate=0.04,
    subsample=0.85,
    colsample_bytree=0.8,
    min_child_weight=3,
    reg_lambda=5.0,
    random_state=0,
)


def log(*a) -> None:
    print(" ".join(map(str, a)), flush=True)


def read_compact(species: str, want: set[str] | None = None) -> list[dict]:
    rows = []
    with (WEB_DIR / f"{species}_compact.jsonl").open() as fh:
        for line in fh:
            r = json.loads(line)
            if want is not None and r["accession"] not in want:
                continue
            rows.append(r)
    return rows


def featurize_species(rows, scale_matrix, coefficients):
    accs, X, y = [], [], []
    for r in rows:
        seq = r["sequence"]
        if len(seq) < 60 or not set(seq) <= _STANDARD_AA:
            continue
        accs.append(r["accession"])
        X.append(featurize_distill(seq, scale_matrix, coefficients))
        y.append(float(r["cg_single"]))
    X = np.nan_to_num(np.array(X), nan=0.0, posinf=0.0, neginf=0.0)
    return np.array(accs), X, np.array(y)


def gate_metrics(pred, web):
    rho = float(st.spearmanr(pred, web).correlation)
    pearson = float(np.corrcoef(pred, web)[0, 1])
    delta = web - pred
    ad = np.abs(delta)
    return {
        "spearman_rho": round(rho, 6),
        "pearson_r": round(pearson, 6),
        "r2": round(float(1 - delta.var() / web.var()), 6),
        "median_abs_delta": round(float(np.median(ad)), 6),
        "mean_abs_delta": round(float(ad.mean()), 6),
        "q90_abs_delta": round(float(np.quantile(ad, 0.90)), 6),
        "n": int(len(web)),
    }


def formula_basis():
    from catgranule._distill import coefficients as coeffs
    from catgranule._distill import scale_matrix as smat
    from catgranule.weights import load_weights

    legacy = load_weights(route="legacy-unaudited")
    return smat(legacy), coeffs(legacy)


def train_save() -> None:
    ts = time.time()
    scale_matrix, coefficients = formula_basis()
    split = json.loads((WEB_DIR / "holdout_human.json").read_text())
    evl, hold = set(split["eval_pool"]), set(split["holdout"])
    rows = read_compact("human", evl | hold)
    accs, X, y = featurize_species(rows, scale_matrix, coefficients)
    is_eval = np.array([a in evl for a in accs])
    Xe, ye = X[is_eval], y[is_eval]
    Xh, yh = X[~is_eval], y[~is_eval]
    log(f"feature dim = {X.shape[1]} (expect 181); eval={Xe.shape[0]} holdout={Xh.shape[0]}")
    assert X.shape[1] == 181, f"unexpected feature dim {X.shape[1]}"
    log(f"featurized in {time.time()-ts:.1f}s")

    model = xgb.XGBRegressor(n_jobs=64, **LOCKED)
    model.fit(Xe, ye)
    pred = model.predict(Xh)
    m = gate_metrics(pred, yh)
    log("HOLDOUT GATE (reproduction):")
    log(json.dumps(m, indent=2))

    rec_rho, rec_med = 0.998961, 0.017653
    ok = abs(m["spearman_rho"] - rec_rho) <= 0.0005 and abs(
        m["median_abs_delta"] - rec_med
    ) <= 0.0005
    log(
        f"reproduce rho {m['spearman_rho']} vs recorded {rec_rho}; "
        f"med|d| {m['median_abs_delta']} vs {rec_med} -> {'PASS' if ok else 'FAIL'}"
    )
    if not ok:
        raise SystemExit(
            "REPRODUCTION GATE FAILED - stopped per D21, no artifacts written"
        )

    model_path = WEIGHTS_DIR / "distilled_xgb_model.ubj"
    model.save_model(str(model_path))
    log(f"saved model -> {model_path} ({model_path.stat().st_size/1e6:.2f} MB)")

    manifest = json.loads((WEB_DIR / "manifest.json").read_text())
    manifest_sha = next(
        i["sha256"] for i in manifest["files"] if i["organism"] == "human"
    )

    from catgranule.weights import load_weights as lw

    legacy = lw(route="legacy-unaudited")
    artifact = {
        "version": "1.0.0",
        "mode": "distilled",
        "model_file": "distilled_xgb_model.ubj",
        "weights": legacy.weights,
        "scales": legacy.scales,
        "normalization": {
            "mean": legacy.normalization_mean,
            "std": legacy.normalization_std,
            "n_proteins": legacy.normalization_n_proteins,
        },
        "window_radius": legacy.window_radius,
        "residue_window_size": legacy.residue_window_size,
        "feature_spec": {
            "dim": 181,
            "map": "window summaries radii 3/7/13/25 (mean/p25/p50/p75/std per 6 "
            "scales=120); g-profile stats half-windows 3/12/25 (11 stats + "
            "frac_pos + mean_pos = 39); AA composition (20); L/logL (2)",
            "filters": "sequence length >= 60, standard residues only",
            "coefficient_basis": "Bolognesi 2016 published coefficients + mmc1 "
            "Table S4 scales (identical copy of the legacy/paper formula "
            "artifacts)",
        },
        "trained_on": {
            "dataset": "PhaSePred web-archive human eval-pool (2022-02-11 snapshot)",
            "n_train": int(is_eval.sum()),
            "human_reviewed.zip_sha256": manifest_sha,
            "compact_source": "human_compact.jsonl <- human_reviewed.json "
            "(data/raw/external/phasepred_web, gitignored)",
            "two_layer_protocol": "config selected on eval-internal split "
            "(train 13000/select 3164, seed 123); holdout scored once",
            "locked_config": LOCKED,
        },
        "provenance": {
            "route": "distilled",
            "source": "XGBoost surrogate fitted to PhaSePred web-archive "
            "catGRANULE single (2022-02-11 snapshot, 16,164 human eval-pool "
            "proteins; 181-dim feature map; locked hyperparameters per "
            "work-log/S4.md breakpoint guide)",
            "date": "2026-09-07",
            "gate": {
                "protocol": "web-archive human holdout (seed 42) single "
                "Spearman rho / median |delta|",
                "rho": m["spearman_rho"],
                "delta_median": m["median_abs_delta"],
                "n_holdout": int(Xh.shape[0]),
                "gate_pass": bool(
                    m["spearman_rho"] >= 0.99 and m["median_abs_delta"] <= 0.01
                ),
                "decision": "accepted with gap by user (D21, 2026-09-07); "
                "original gate delta<=0.01 superseded",
                "decision_ref": "docs/migration/PLAN.md D21",
            },
        },
    }
    art_path = WEIGHTS_DIR / "catgranule_distilled_weights.json"
    art_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    log(f"saved artifact -> {art_path}")

    report = {
        "reproduced_numbers": m,
        "recorded_numbers": {
            "spearman_rho": rec_rho,
            "median_abs_delta": rec_med,
        },
        "model_file": str(model_path),
        "model_bytes": model_path.stat().st_size,
        "artifact_file": str(art_path),
        "training": {
            "n_train": int(is_eval.sum()),
            "n_holdout": int(Xh.shape[0]),
            "locked_config": LOCKED,
        },
    }
    (WEB_DIR / "s4b_reproduce_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    log("wrote s4b_reproduce_report.json to data area")


def paper_artifact() -> None:
    cand = json.loads(
        (WEB_DIR / "candidates/catgranule_paper_route_candidate.json").read_text()
    )
    norm = json.loads((WEB_DIR / "paper_route_yeast_norm.json").read_text())
    cand["version"] = "1.0.0"
    cand["normalization"] = {
        "mean": norm["mean"],
        "std": norm["std"],
        "n_proteins": norm["n_proteins"],
    }
    cand["provenance"] = {
        "route": "paper",
        "source": "Bolognesi 2016 (mmc7 Eq1-2 published coefficients "
        "aR=(0.48,7.24) aD=(0.26,11.54) aP=(1.98,1.42) aL=0.25; mmc1 Table S4 "
        "scales); yeast proteome Z-normalization recomputed "
        f"2026-09-06 (data/interim/uniprot_cache.jsonl, n={norm['n_proteins']}); "
        "deterministic closed-form",
        "date": "2026-09-07",
        "gate": {
            "rho": 0.974712,
            "delta_median": 0.186215,
            "holdout": "web-human (seed 42, n=4070)",
            "pass": False,
            "position": "alternate artifact; original S4 main route (below "
            "gate, see CATGRANULE_VALIDATION.md)",
        },
    }
    out = WEIGHTS_DIR / "catgranule_paper_weights.json"
    out.write_text(json.dumps(cand, indent=2), encoding="utf-8")
    log(f"saved paper artifact -> {out}")


def cross_species() -> None:
    from catgranule.surrogate import load_surrogate
    from catgranule.weights import load_weights

    weights = load_weights(route="distilled")
    eng = load_surrogate(weights)
    scale_matrix, coefficients = formula_basis()
    out = {}
    for species in ("mouse", "yeast"):
        rows = read_compact(species)
        accs, X, y = featurize_species(rows, scale_matrix, coefficients)
        pred = eng.predict(X)
        m = gate_metrics(pred, y)
        tot = len(rows)
        out[species] = {
            **m,
            "archived_total": tot,
            "scored": int(len(pred)),
            "coverage": round(float(len(pred)) / tot, 6),
        }
        log(
            f"CROSS-SPECIES {species}: scored {len(pred)}/{tot} ({len(pred)/tot:.3%}) "
            f"rho={m['spearman_rho']} med|d|={m['median_abs_delta']}"
        )
    (WEB_DIR / "s4b_cross_species.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    log("wrote s4b_cross_species.json to data area")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-save", action="store_true")
    ap.add_argument("--paper-artifact", action="store_true")
    ap.add_argument("--cross-species", action="store_true")
    args = ap.parse_args()
    if args.train_save:
        train_save()
    if args.paper_artifact:
        paper_artifact()
    if args.cross_species:
        cross_species()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())