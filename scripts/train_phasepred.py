#!/usr/bin/env python
"""One-command retraining of the four PhaSePred models from the committed tables.

Reads data/processed/*.tsv (E1), applies the paper protocol (Chen et al. 2022
PNAS, RETRAIN_PROTOCOL.md), and writes 40 XGBClassifier artifacts plus their
manifests, metrics.json, leakage_report.json and the sealed accession lists.

Reproduction is byte-exact and depends on xgboost 3.2.x; see PLAN-E2E §2.9.

Usage:
    python scripts/train_phasepred.py                          # all four tasks
    python scripts/train_phasepred.py --task SaPS --task PdPS  # subset
    python scripts/train_phasepred.py --smoke 200 --output-dir /tmp/smoke
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
import xgboost
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from phasepred.data import data_path, models_root  # noqa: E402
from phasepred.training import (  # noqa: E402
    N_FOLDS,
    N_MODELS,
    NEG_RATIO,
    SEED,
    TASKS,
    check_leakage,
    ensemble_score,
    feature_columns,
    run_cross_validation,
    split_pos_neg,
    task_frame,
    train_final_ensemble,
)
from phasepred.training_data import FILENAMES, sha256_file, table_path  # noqa: E402

FEATURE_DEFINITIONS_VERSION = "v2022"
PROTOCOL = "paper-default"


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task",
        action="append",
        default=None,
        metavar="TASK",
        help="task to train, repeatable or comma-separated "
        f"(one of {', '.join(TASKS)}; default: all)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="directory for model artifacts (default: the bundled models root)",
    )
    parser.add_argument(
        "--metrics-dir",
        type=Path,
        default=None,
        help="directory for metrics.json/manifest.json (default: the product dir)",
    )
    parser.add_argument(
        "--smoke",
        type=int,
        default=None,
        metavar="N",
        help="use only the first N negatives per task; implies --skip-cv and "
        "requires explicit --output-dir and --metrics-dir outside the published paths",
    )
    parser.add_argument(
        "--skip-cv",
        action="store_true",
        help="skip the 5-fold cross-validation",
    )
    return parser.parse_args(argv)


def _select_tasks(raw: list[str] | None) -> list[str]:
    tasks: list[str] = []
    for value in raw or []:
        tasks.extend(part.strip() for part in value.split(",") if part.strip())
    return tasks or list(TASKS)


def _is_within(path: Path, root: Path) -> bool:
    path = Path(path).resolve()
    root = Path(root).resolve()
    return path == root or root in path.parents


def resolve_output_dirs(
    output_dir: Path | None,
    metrics_dir: Path | None,
    smoke: int | None,
) -> tuple[Path, Path]:
    """Resolve where models and metrics go, refusing unsafe smoke targets.

    Smoke mode truncates the negative pool, so its metrics and accession lists
    are not the published ones. It may not write either default location.
    """
    published_models = models_root()
    published_metrics = data_path("products", "A_paper_split_recomputed")
    resolved_models = Path(output_dir) if output_dir is not None else published_models
    resolved_metrics = Path(metrics_dir) if metrics_dir is not None else published_metrics
    if smoke is None:
        return resolved_models, resolved_metrics

    refusals: list[str] = []
    published_roots = (
        (published_models, "models"),
        (published_metrics, "product"),
    )
    for target, flag, label in (
        (resolved_models, output_dir, "--output-dir"),
        (resolved_metrics, metrics_dir, "--metrics-dir"),
    ):
        if flag is None:
            refusals.append(
                f"--smoke needs an explicit {label}; refusing to overwrite "
                f"{target} (smoke results are truncated)"
            )
            continue
        for root, kind in published_roots:
            if _is_within(target, root):
                refusals.append(
                    f"--smoke refuses to overwrite the published {kind} directory "
                    f"{target} (smoke results are truncated)"
                )
                break
    if refusals:
        for refusal in refusals:
            print(f"error: {refusal}", file=sys.stderr)
        raise SystemExit(2)
    return resolved_models, resolved_metrics


def _run_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT), capture_output=True, text=True
        )
        return out.stdout.strip() if out.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def _table_record(scope: str, split: str) -> dict[str, str]:
    return {
        "path": str(Path("data") / "processed" / FILENAMES[(scope, split)]),
        "sha256": sha256_file(table_path(scope, split)),
    }


def _input_tables() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for scope, split in FILENAMES:
        records[f"{scope}_{split}"] = _table_record(scope, split)
    return records


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def _print_environment(commit: str) -> None:
    print(f"train_phasepred — commit {commit}")
    print(
        f"  python {sys.version.split()[0]}; xgboost {xgboost.__version__}; "
        f"scikit-learn {sklearn.__version__}; numpy {np.__version__}; pandas {pd.__version__}"
    )
    for scope, split in FILENAMES:
        record = _table_record(scope, split)
        print(f"  {record['path']}  {record['sha256']}")


def _train_task(
    task: str,
    output_dir: Path,
    smoke: int | None,
    skip_cv: bool,
    trained_at: str,
) -> tuple[dict[str, Any], dict[str, Any], str, str]:
    scope = "human" if task.startswith("h") else "base"
    train_frame = task_frame(task, "train")
    test_frame = task_frame(task, "test")

    overlap = check_leakage(train_frame, test_frame)
    if overlap:
        raise RuntimeError(f"Leakage: {overlap} accessions in both train and test for {task}")

    cols = feature_columns(task)
    X_pos, X_neg_pool = split_pos_neg(train_frame, cols)
    if smoke is not None:
        X_neg_pool = X_neg_pool.head(smoke).reset_index(drop=True)
    n_neg_train = int(len(X_neg_pool))

    print(f"\n=== Training {task} ===")
    print(f"  pos={len(X_pos)} neg_pool={len(X_neg_pool)} test_rows={len(test_frame)}")

    cv_auc = float("nan")
    fold_aucs: list[float] = []
    cv_ens = float("nan")
    fold_ens: list[float] = []
    if not skip_cv:
        cv_auc, fold_aucs, cv_ens, fold_ens = run_cross_validation(X_pos, X_neg_pool, cols, SEED)
        print(
            f"  CV AUC per-model mean (paper-comparable): {cv_auc:.4f}; "
            f"ensemble: {cv_ens:.4f}"
        )

    ensemble = train_final_ensemble(X_pos, X_neg_pool, cols, SEED)
    models_out = output_dir / task
    models_out.mkdir(parents=True, exist_ok=True)
    for i, model in enumerate(ensemble):
        joblib.dump(model, models_out / f"8f_model_{i}.joblib")

    X_test = test_frame[cols].astype(float)
    y_test = test_frame["label"].astype(int).to_numpy()
    test_auc = (
        float(roc_auc_score(y_test, ensemble_score(ensemble, X_test)))
        if len(set(y_test.tolist())) > 1
        else float("nan")
    )
    print(f"  TEST AUC (S3): {test_auc:.4f}  (n_test={len(X_test)})")

    manifest = {
        "task": task,
        "feature_columns": cols,
        "feature_definitions_version": FEATURE_DEFINITIONS_VERSION,
        "protocol": PROTOCOL,
        "n_models": N_MODELS,
        "n_neg_ratio": NEG_RATIO,
        "n_folds": N_FOLDS,
        "seed": SEED,
        "train_sheet": "h" + task[1:] if task.startswith("h") else task,
        "trained_at": trained_at,
        "xgboost_version": xgboost.__version__,
        "auc_cv": cv_auc,
        "auc_cv_ensemble": cv_ens,
        "auc_test_s3": test_auc,
        "input_tables": {
            "train": _table_record(scope, "train"),
            "test": _table_record(scope, "test"),
        },
    }
    _write_json(models_out / "manifest.json", manifest)

    metrics = {
        "train_rows": int(len(train_frame)),
        "test_rows": int(len(test_frame)),
        "n_pos_train": int(len(X_pos)),
        "n_neg_train": n_neg_train,
        "cv_auc": cv_auc,
        "cv_auc_ensemble": cv_ens,
        "cv_fold_aucs": [round(a, 6) for a in fold_aucs],
        "cv_fold_aucs_ensemble": [round(a, 6) for a in fold_ens],
        "test_auc_S3": test_auc,
        "n_neg_per_model": min(n_neg_train, int(len(X_pos) * NEG_RATIO)),
    }
    train_accs = "\n".join(sorted(set(train_frame["UniprotEntry"].astype(str)))) + "\n"
    test_accs = "\n".join(sorted(set(test_frame["UniprotEntry"].astype(str)))) + "\n"
    return metrics, {"train_test_overlap": int(overlap)}, train_accs, test_accs


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    tasks = _select_tasks(args.task)
    unknown = [task for task in tasks if task not in TASKS]
    if unknown:
        print(f"error: unknown task(s) {unknown}; expected one of {list(TASKS)}", file=sys.stderr)
        return 2

    output_dir, metrics_dir = resolve_output_dirs(args.output_dir, args.metrics_dir, args.smoke)

    skip_cv = args.skip_cv
    if args.smoke is not None:
        skip_cv = True
        print(f"--smoke {args.smoke}: cross-validation skipped (smoke negatives are truncated)")

    try:
        trained_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        commit = _run_commit()
        _print_environment(commit)
        input_tables = _input_tables()
        metrics_dir.mkdir(parents=True, exist_ok=True)

        metrics: dict[str, Any] = {}
        leakage: dict[str, Any] = {}
        for task in tasks:
            task_metrics, task_leakage, train_accs, test_accs = _train_task(
                task, output_dir, args.smoke, skip_cv, trained_at
            )
            metrics[task] = task_metrics
            leakage[task] = task_leakage
            (metrics_dir / f"train_accessions_{task}.tsv").write_text(train_accs)
            (metrics_dir / f"test_accessions_{task}.tsv").write_text(test_accs)

        summary = {
            "feature_definitions_version": FEATURE_DEFINITIONS_VERSION,
            "protocol": PROTOCOL,
            "seed": SEED,
            "n_models": N_MODELS,
            "n_neg_ratio": NEG_RATIO,
            "n_folds": N_FOLDS,
            "trained_at": trained_at,
            "xgboost_version": xgboost.__version__,
            "commit": commit,
            "input_tables": input_tables,
        }
        _write_json(metrics_dir / "metrics.json", metrics)
        _write_json(metrics_dir / "leakage_report.json", leakage)
        _write_json(metrics_dir / "manifest.json", summary)

        print(f"\nWrote {len(tasks) * N_MODELS} models under {output_dir}")
        print(f"Wrote metrics.json, manifest.json, leakage_report.json under {metrics_dir}")
    except Exception as exc:  # noqa: BLE001
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
