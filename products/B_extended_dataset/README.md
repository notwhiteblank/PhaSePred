# Product B — Extended dataset + tuned XGBoost (default)

This is the default production model. It departs from Product A in two
ways:

1. **Training labels come from the current databases**: PhaSepDB 3.0
   (positive labels for SaPS / PdPS), the paper's 10-organism proteome
   pool (negatives), and LLPSDB v2 confirmed negatives.
2. **Imbalance handled with `scale_pos_weight`, not 2:1 subsampling**.
   The negative pool is now ~60,000; subsampling discards too much signal.

Hyperparameters are again Optuna-tuned per task on the training partition
only. Optuna 5-fold inner CV; final evaluation on a sealed 20% test
partition.

## Files

| File | Meaning |
|---|---|
| `tune.py` | Re-run Optuna (100 trials × 4 tasks) on the training partition. |
| `train.py` | Train the 10-model ensemble using `tuned_params.json`. |
| `tuned_params.json` | Per-task best params. |
| `train_accessions.tsv` / `test_accessions.tsv` | The sealed split. |
| `leakage_report.json` | Critical: PhaSepDB 3.0 overlaps PhaSePro, so external AUC on PhaSePro must filter overlapping accessions. This file records what was removed. |
| `metrics.json` | CV AUC, held-out test AUC, external AUCs on non-overlapping subsets. |
| `models/<task>/8f_model_{0..9}.joblib` | Trained pipelines. |

## Reproduce from scratch

```bash
uv run python products/B_extended_dataset/tune.py
uv run python products/B_extended_dataset/train.py
```

## Use this product

```bash
uv run phasepred predict --fasta my.fasta --mode SaPS
# Product B is the default; no --product flag needed.
```
