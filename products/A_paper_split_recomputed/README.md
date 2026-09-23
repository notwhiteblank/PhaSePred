# Product A — paper-protocol retrain (v2022, Dataset S2/S3 sheet features)

Retrained per the faithful paper protocol,
XGBoost with **library-default hyperparameters**, trained on the Dataset S2/S3
sheet feature columns (not on locally recomputed tool outputs), with the
paper's 5-fold × 5-round × 10 negative-set (2:1) ensemble scheme.

## Files

| File | Meaning |
|---|---|
| `scripts/train_phasepred.py` | Canonical retrain script (paper protocol, SEED 42). `products/A_paper_split_recomputed/train.py` is a thin wrapper around it. |
| `scripts/eval_paper_layer.py` | Paper-layer six-row CV table vs mineru_md:67 (`--mode paper` for sheet features, `--mode control --features <csv>` for the locally-recomputed control arm). |
| `src/phasepred/data/models/<task>/8f_model_{0..9}.joblib` | 10 raw `XGBClassifier` models per task (6/8 features all-species, 10 human). The trained models are shipped as package data under `src/`; this directory holds the run provenance and splits, not a `models/` tree. |
| `src/phasepred/data/models/<task>/manifest.json` | Per-mode manifest binding `feature_definitions_version: v2022` (stale-warning contract). |
| `manifest.json` | Run provenance (inputs sha256, versions, commit). |
| `metrics.json` | CV AUC (fold-means), S3 test AUC, auxiliary external AUCs, `feature_definitions_version`. |
| `leakage_report.json` | S2-train / S3-test accession overlap (must be 0). |
| `train_accessions_*.tsv`, `test_accessions_*.tsv` | Accession lists per task (S2 / S3). |
| `tuned_params.json`, `tune.py` | Historical Optuna control (pre-v2022, kept as reference only; not used by the v2022 artifacts). |

## Reproduce from scratch

```bash
uv run python scripts/train_phasepred.py
uv run python scripts/eval_paper_layer.py --mode paper --output paper_layer_cv_results.json
```

## Use this product

```bash
uv run phasepred predict --fasta my.fasta --mode SaPS --output scores.csv
```
