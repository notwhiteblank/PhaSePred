# Product A — Paper split + recomputed features + tuned XGBoost

This product reproduces the paper's training protocol but with two
deliberate departures from the original PhaSePred:

1. **Features are recomputed by our own tool runs**, not taken from the
   paper supplement. This makes the model end-to-end reproducible from
   FASTA without depending on the paper's intermediate spreadsheets.
2. **XGBoost hyperparameters are Optuna-tuned** on the training partition
   (S2), then frozen. The original paper used default XGBoost params.

Everything else follows the paper:

- Train / test split = paper's S2 / S3
- 10-model ensemble per task, each model fitted on a 2:1 negative subsample
  of S2
- Final prediction = mean over the 10 models

## Files

| File | Meaning |
|---|---|
| `tune.py` | Re-run Optuna (100 trials × 4 tasks). Writes `tuned_params.json`. |
| `train.py` | Train the 10-model ensemble using `tuned_params.json`. Writes `models/`, `metrics.json`, `train_accessions.tsv`, `test_accessions.tsv`, `leakage_report.json`. |
| `tuned_params.json` | Per-task best params (n_estimators, max_depth, learning_rate, …). |
| `train_accessions.tsv` | UniProt accessions used for training (S2). |
| `test_accessions.tsv` | UniProt accessions held out for evaluation (S3). |
| `leakage_report.json` | Overlap counts between training and PhaSePro / LLPSDB2. |
| `metrics.json` | CV AUC, S3 test AUC, external AUCs (on non-overlapping subset). |
| `models/<task>/8f_model_{0..9}.joblib` | Trained sklearn pipelines. |

## Reproduce from scratch

```bash
uv run python products/A_paper_split_recomputed/tune.py
uv run python products/A_paper_split_recomputed/train.py
```

## Use this product

```bash
uv run phasepred predict --fasta my.fasta --mode SaPS --product A
```
