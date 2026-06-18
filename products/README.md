# Products

Two shippable, independent PhaSePred models.

| Product | Data | Features | Model | Default in `predict`? |
|---|---|---|---|---|
| **A** `A_paper_split_recomputed/` | Paper S2/S3 split | Recomputed by our own tool runs (`data/interim/recomputed_features_full.csv`) | XGBoost 10-model ensemble, Optuna-tuned, 2:1 negative subsampling (paper protocol) | no |
| **B** `B_extended_dataset/` | PhaSepDB 3.0 + paper's 10-organism negatives, 80/20 stratified split | Recomputed by our own tool runs | XGBoost 10-model ensemble, Optuna-tuned, `scale_pos_weight` (all negatives) | **yes** |

Each product directory contains:

```
<product>/
  README.md                       How to use this product
  tune.py                         Run Optuna tuning from scratch
  train.py                        Train the 10-model ensemble using tuned params
  tuned_params.json               Best XGBoost hyperparameters per task
  train_accessions_<task>.tsv     Sealed training-set UniProt accessions per task
  test_accessions_<task>.tsv      Sealed held-out-set UniProt accessions per task
  leakage_report.json             Overlap analysis between training labels and
                                  external validation sets (PhaSePro, LLPSDB2)
  metrics.json                    CV AUC, held-out test AUC, external AUC (on the
                                  non-overlapping subset only)
  models/<task>/8f_model_*.joblib 10 trained XGBoost pipelines per task
```

`uv run phasepred predict` defaults to Product B. Use `--product A` to
switch.

## Anti-leakage protocol

Both products follow the same rules:

1. Train / test split is computed once, deterministically (seeded), and
   the test set is sealed before tuning starts.
2. Optuna does 5-fold stratified CV on **training data only**; the test
   set is never visible to the tuner.
3. External validation sets (PhaSePro, LLPSDB2) are filtered to the
   subset of accessions that do **not** appear in training. AUCs are
   reported on that non-overlapping subset; the overlap counts and the
   exact accessions removed are recorded in `leakage_report.json`.
4. The training pipeline aborts if the same accession ends up in both
   train and test files.
