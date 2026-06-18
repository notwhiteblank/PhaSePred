# PhaSePred

> Predict the phase-separation propensity of a protein from its amino-acid
> sequence. Maintained by the **Tingting Li Lab** at Peking University.
> This repository is the updated successor to our group's earlier work
> [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119) —
> retrained on the expanded current databases with Optuna-tuned XGBoost,
> packaged as an end-to-end CLI.

[中文](README.md)

## Background

**Phase separation** is the spontaneous condensation of proteins into
liquid droplet-like compartments — picture oil droplets in water. These
membraneless organelles (stress granules, the nucleolus, Cajal bodies,
…) organize biochemistry without a lipid bilayer. Misregulated phase
separation is implicated in ALS, Alzheimer's, and several cancers.

PhaSePred scores a protein's tendency to participate in phase separation
under one of four task definitions:

| Mode | Meaning | Species scope |
|------|------|---------|
| **SaPS** | Self-assembling phase separation — the protein condenses on its own | any (8 features) |
| **PdPS** | Partner-dependent phase separation — needs interaction partners | any (8 features) |
| **hSaPS** | Human-only SaPS | human (10 features) |
| **hPdPS** | Human-only PdPS | human (10 features) |

The output is a score in [0, 1]: higher means more likely to belong to
that mode. It is **not** a calibrated probability, only a relative
ranking within the mode.

## Two independent products

The repository ships **two models**, both shippable as standalone
software:

| Product | Training data | Split | Class imbalance | Default? |
|---|---|---|---|---|
| **A** | Paper S2 / S3 + our recomputed features | Paper's original split | 2:1 negative subsample (paper protocol) | no |
| **B** | PhaSepDB 3.0 + 10-organism background negatives + LLPSDB v2 + recomputed features | Stratified 80/20, seed=42 | `scale_pos_weight` (all negatives) | **yes** |

Both use Optuna for 100-trial × 5-fold CV tuning on the training
partition. **The test partition is sealed**: tuning never sees it.
`predict` defaults to Product B; `--product A` selects A. See
[products/README.md](products/README.md) for the full spec.

### Key numbers

| Task | A: CV → S3 test | B: CV → 20% test | B: PhaSePro external |
|------|------|------|------|
| SaPS  | 0.888 → **0.794** | 0.840 → **0.838** | 0.782 |
| PdPS  | 0.759 → **0.699** | 0.764 → **0.737** | 0.859 |
| hSaPS | 0.933 → **0.843** | 0.855 → **0.816** | 0.839 |
| hPdPS | 0.851 → **0.785** | 0.823 → **0.859** | 0.959 |

> External validation sets (PhaSePro, LLPSDB2) are **filtered to drop any
> accession present in the training set**. Counts of the removed entries
> and the surviving subset sizes are written to
> `products/<product>/leakage_report.json`. The LLPSDB2 column hovers
> around 0.5 for both products — its "confirmed negatives" sit near the
> decision boundary, which limits binary discrimination. This is a label
> property, not a model defect; recorded for honesty.

## Quick start

```bash
# 1. Clone and bring up the Python env
git clone https://github.com/notwhiteblank/PhaSePred.git && cd PhaSePred
uv sync

# 2. Install the external feature tools (next section)
bash tools/install/install_all.sh    # or follow the per-tool steps below

# 3. Verify everything is found
uv run phasepred check-tools

# 4. Predict
uv run phasepred predict --fasta my.fasta --mode SaPS --output scores.csv
```

A batch takes ~8 seconds. The output CSV includes every computed feature
column so you can sanity-check the inputs.

## Installing the external tools (required)

PhaSePred invokes **6 external feature tools** (catGRANULE is reimplemented
as a built-in Python module — no install needed). Every tool is located
through the path resolver — order is `env var → vendored in repo → PATH`.
After installing, always run `phasepred check-tools` to verify.

| Tool | License | Install | Env override |
|---|---|---|---|
| **SEG** | NCBI public-domain convention | Vendored; `cd tools/per-tool/SEG && make` | `PHASEPRED_SEG_BIN` |
| **PScore** | CC-BY 4.0 | Vendored; no install needed | `PHASEPRED_PSCORE_DIR` |
| **PLAAC** | MIT | Vendored prebuilt jar; needs Java 17+ | `PHASEPRED_PLAAC_JAR` |
| **ESpritz** | Academic license, no redistribution | `bash tools/install/install_espritz.sh` (prompts you to accept the upstream license) | `PHASEPRED_ESPRITZ_DIR` |
| **IUPred3** | Academic license, no redistribution | `bash tools/install/install_iupred3.sh` (same) | `PHASEPRED_IUPRED3_DIR` |
| **DeepCoil** | No upstream license declared | `bash tools/install/install_deepcoil_env.sh` (Python 3.8 conda env) | `PHASEPRED_DEEPCOIL_ENV` |

The vendored tools just work after `git clone`. Tools forbidden from
redistribution must be downloaded by the user. The full license audit is
at [docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md).

### How the path resolver works

There is **no `Tools/<X>/<bin>` hardcoded anywhere** in the source.
`src/phasepred/tool_paths.py` exposes `find_seg() / find_pscore_dir() /
find_plaac_wrapper() / …`, each running a 4-step cascade:

1. **Environment variable**, e.g. `PHASEPRED_SEG_BIN=/opt/bin/seg`
2. **Vendored in repo**: `tools/per-tool/SEG/seg`
3. **PATH**: `shutil.which("seg")`
4. **Error**: raises `PhaSePredToolNotFound` with the install hint string

If you already have a tool installed system-wide (conda, apt, manual
build), **don't reinstall and don't copy it into the repo** — just point
the env var at it.

### The check-tools workflow

```bash
uv run phasepred check-tools
```

Example output:

```
TOOL        STATUS  REQUIRED  SOURCE                       PATH
SEG         OK      yes       vendored                     /repo/tools/per-tool/SEG/seg
PScore      OK      yes       vendored                     /repo/tools/per-tool/PScore/SourceCodeS2
PLAAC       OK      yes       vendored                     /repo/tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar
IUPred3     OK      no        env:PHASEPRED_IUPRED3_DIR    /opt/iupred3
ESpritz     OK      yes       env:PHASEPRED_ESPRITZ_DIR    /opt/espritz
DeepCoil    OK      yes       vendored (conda)             .external_envs/deepcoil
catGRANULE  OK      yes       built-in                     src/phasepred/catgranule_v1.py

7 OK  ·  0 MISSING  ·  0 required-but-missing
```

If anything is `MISSING` and marked `required: yes`, `phasepred predict`
exits **before** reading the FASTA with the same status table and a
per-tool install hint — no half-completed pipeline.

### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `SEG executable not found` | SEG wasn't built | `cd tools/per-tool/SEG && make` |
| `PLAAC: java not found` | No Java runtime | `sudo apt install openjdk-17-jre-headless`, or macOS `brew install openjdk` |
| `ESpritz: perl not found` | No Perl | Usually preinstalled on Linux; macOS `brew install perl` |
| `DeepCoil env not found` | Conda env not created | `bash tools/install/install_deepcoil_env.sh` |
| `IUPred3 not found` after download | Wrong extraction path | Put `iupred3/` at `tools/per-tool/IUPred3/iupred3/` or set `PHASEPRED_IUPRED3_DIR` |

## Optional: extra data for human-mode predictions

`SaPS` and `PdPS` use 8 features that are all computed from sequence — no
extra files needed.

`hSaPS` and `hPdPS` add 2 features that come from upstream **lookup
tables**, not from inference tools: `Phos freq` (phosphorylation-site
density) and `DeepPhase` (microscopy-based phase-separation score).
Neither has a public inference model — PhaSePred **cannot recompute
them** from sequence. The two files must be downloaded by the user:

| Feature | Source | License | Expected path |
|---|---|---|---|
| Phos freq | [PhosphoSitePlus](https://www.phosphosite.org/) `Phosphorylation_site_dataset.gz` | registration required, no redistribution | `data/raw/external/phosphositeplus/Phosphorylation_site_dataset.gz` |
| DeepPhase | [DeepPhase paper supplement](https://github.com/PEILab/DeepPhase) `tableS3.xlsx` | academic, see upstream | `data/raw/external/deepphase/extracted/tableS3.xlsx` |

**What happens if the files are missing?** `phasepred predict --mode
hSaPS|hPdPS` **still runs**, but emits explicit
`PhaSePredMissingDataWarning` entries on stderr telling you exactly
which proteins were median-imputed (and therefore scored on partial
features):

```
PhaSePredMissingDataWarning: DeepPhase has no entry for 3/5 requested
protein(s): P12345, Q67890, A1B2C3. Median imputation will be used, so
hSaPS/hPdPS scores for these proteins should be treated as approximate.
```

If you only need SaPS/PdPS predictions you can ignore both files.

## Usage

### Predict (FASTA → scores)

```bash
# 8-feature models (any species)
uv run phasepred predict --fasta proteins.fasta --mode SaPS  --output scores.csv
uv run phasepred predict --fasta proteins.fasta --mode PdPS  --output scores.csv

# 10-feature models (human only)
uv run phasepred predict --fasta proteins.fasta --mode hSaPS --output scores.csv
uv run phasepred predict --fasta proteins.fasta --mode hPdPS --output scores.csv

# Switch product
uv run phasepred predict --fasta proteins.fasta --mode SaPS --product A --output scores.csv

# Predict from UniProt IDs (three input forms)
uv run phasepred predict --ids "P35637,Q9Y2W1" --mode SaPS --output scores.csv      # inline
uv run phasepred predict --ids my_ids.txt --mode SaPS --output scores.csv           # file (one ID per line)
cat my_ids.txt | uv run phasepred predict --ids - --mode SaPS --output scores.csv   # stdin

# FASTA + IDs can be combined (deduplicated after merge)
uv run phasepred predict --fasta my.fasta --ids "P35637" --mode SaPS --output scores.csv
```

ID input calls `https://rest.uniprot.org` to fetch sequences; results
are cached at `data/interim/uniprot_cache.jsonl`, so repeat predictions
on the same IDs do not hit the network. Override the cache path with
`--cache PATH`.

### Reproduce both products

```bash
# 1. Prepare features (only needed if data/interim/recomputed_features_full.csv is missing).
#    Requires every external tool to be installed.
uv run phasepred features-recomputed \
    --sequences data/interim/sequences.csv \
    --output data/interim/recomputed_features_full.csv \
    --espritz-cache data/interim/espritz_idr_cache.jsonl \
    --catgranule-source v1

# 2. Tune + train Product A
uv run python products/A_paper_split_recomputed/tune.py
uv run python products/A_paper_split_recomputed/train.py

# 3. Tune + train Product B
uv run python products/B_extended_dataset/tune.py
uv run python products/B_extended_dataset/train.py
```

Tuning takes 30-60 minutes per product; training 1-2 minutes.

## Repository layout

```
src/phasepred/             Shared library + CLI
  cli.py                   predict / check-tools / features-recomputed / ...
  predictor.py             FASTA → features → score pipeline; --product A|B
  tool_paths.py            External-tool resolver (env / vendored / PATH)
  features.py              Native features (Hydropathy, FCR) and schema
  catgranule_v1.py         Paper-formula catGRANULE reimplementation
  paper.py                 Used by Product A: paper S2/S3 supplement parser
  updated_data.py          Used by Product B: PhaSepDB3 + LLPSDB2 merger
  ...

products/                  Two independent product deliverables
  A_paper_split_recomputed/  Paper split + recomputed features + tuned XGB
    tune.py / train.py
    tuned_params.json metrics.json leakage_report.json
    train_accessions_*.tsv test_accessions_*.tsv   (sealed partitions)
    models/<task>/8f_model_{0..9}.joblib            (10-model ensembles)
  B_extended_dataset/        Extended dataset + tuned XGB (default)
    (same shape as above)

tools/                     External tool layer
  README.md                            Resolver mechanics
  THIRD_PARTY_LICENSES.md              License audit
  install/                             install_<tool>.sh scripts
  wrappers/                            run_<tool>.sh shell wrappers
  per-tool/<TOOL>/                     Per-tool README and (where licensed) vendored source

docs/
  THIRD_PARTY_LICENSES.md
  DATA_SOURCES.md
  ENVIRONMENT.md
  REPO_LAYOUT.md
  STORAGE_LAYOUT.md

tests/                     pytest suite
```

## The exploration narrative and the numbers behind it

The project began as an internal rebuild of the Tingting Li Lab's 2022
PhaSePred codebase (Chen et al., *PNAS*) and evolved, through database
updates and methodology polish, into the two products shipped in this
release. Each step is recorded here with **the concrete numbers it
produced** — readers don't need any additional documents to understand
what happened.

### Step 1a · LogReg sanity check

LogisticRegression on three simple features (Hydropathy, FCR, IDR) using
the 2022 S2/S3 split. Goal: validate the data-loading and evaluation
plumbing. AUCs landed well below the 2022 published values
(~0.7-range), which was expected — too few features. Pipeline confirmed
working.

### Step 1b · XGBoost Optuna exploration

100-trial × 3-fold CV Optuna search over 9 XGBoost hyperparameters
(n_estimators 100-600, max_depth 3-8, log-uniform learning_rate
0.01-0.3, etc.) on the 2022 paper-supplied feature values. Exploration
only; no final model produced.

### Step 1c · End-to-end rebuild of the 2022 protocol

Paper-supplied features + default XGBoost + 2:1 negative subsampling +
10-model ensemble — strict rebuild of the 2022 Table 1 to confirm the
new code stack faithfully reproduces our prior published numbers:

| Task | 2022 published | This repo's rebuild (S3 test AUC) |
|---|---|---|
| SaPS  | 0.86 | 0.79 |
| PdPS  | 0.66 | 0.74 |
| hSaPS | 0.85 | 0.83 |
| hPdPS | 0.78 | 0.82 |

All four tasks within ±0.05 of the 2022 published numbers; the residual
variance is from the 2:1 subsampling. This confirms the new code stack
faithfully reproduces the lab's original methodology.

### Step 2 · External tool integration

Installed and validated SEG (NCBI public-domain), ESpritz (academic
license), IUPred3 (academic license), PScore (CC-BY 4.0), PLAAC (MIT),
DeepCoil (no upstream license), and catGRANULE so predictions could run
straight from FASTA without the 2022 paper supplements. CD-HIT was
installed during this phase as a sequence-redundancy filter but never
integrated into either product; removed.

### Step 3 · AlphaFold pLDDT feature ablation (NEGATIVE result)

**Hypothesis**: AlphaFold's per-residue pLDDT mean carries
structured/disordered information missing from the existing features.
Adding it as a 9th feature should improve AUC.

**Method**: Pulled mean pLDDT for 60,315 proteins via the AlphaFold DB
API; retrained the 2022 protocol ensemble with 9 features; compared.

**Result (negative)**:

| Task | 8 features AUC | 9 features (with pLDDT) AUC | Δ |
|---|---|---|---|
| SaPS  | 0.794 | 0.794 |  0.000 |
| PdPS  | 0.699 | 0.701 | +0.002 |
| hSaPS | 0.843 | 0.838 | -0.005 |
| hPdPS | 0.785 | 0.781 | -0.004 |

Feature importance ranked pLDDT 8 of 9, with substantial correlation to
PScore (r ≈ 0.45), LCR (r ≈ 0.50), and PLAAC (r ≈ 0.42) — pLDDT's
"structured-vs-disordered" signal was already captured. **Decision: drop
the feature.** Production uses 8 base features (10 for human-only models
that add Phos freq + DeepPhase). A useful negative result: future work
should look elsewhere.

### Step 4 · Dataset-expansion first pass

After PhaSepDB 3.0 was released, redid the 2022 split protocol
(PhaSepDB3 positives + the 10-organism background + LLPSDB v2 confirmed
negatives) and recomputed every feature. **Issue**: this iteration used
default XGBoost params and **did not** filter accession overlap on
external validation (PhaSePro / LLPSDB2). The 5-fold CV AUCs looked
excellent (SaPS 0.866, PdPS 0.832, hSaPS 0.977, hPdPS 0.973) but the
optimistic external numbers exposed the missing overlap filter.
Superseded by Product B, which adds Optuna tuning and proper accession
filtering.

### Auxiliary · catGRANULE v1 reconstruction

catGRANULE v1's upstream code was unavailable, only the formula from
Bolognesi et al. 2016. Reimplemented from scratch: parsed the paper's
supplementary Table S1 to derive Z-normalization constants, coded the
granule-propensity score by formula. Validated against the catGRANULE
column in our 2022 PhaSePred supplement S2: **RMSE ≈ 0.24**, feature-
level AUC at parity. The pipeline prefers 2022 supplement values when
available and falls back to this reimplementation
(`src/phasepred/catgranule_v1.py`) otherwise.

### Step 5 · IDR feature sensitivity (ESpritz vs IUPred3)

The 2022 paper specifies ESpritz DisProt at 5% FPR for the IDR signal;
IUPred3 is a more common drop-in. Held everything else constant and
swapped the IDR source: ESpritz beat IUPred3 by 0.01-0.03 AUC across
all four tasks, and matched the 2022 definition. **Kept ESpritz** as
the default; IUPred3 remains available in `src/phasepred/iupred.py` as
an alternative.

### This release · Products A and B (final deliverables)

Both products use:
- **Optuna 100 trials × 5-fold stratified CV** on the training partition
  only (per-task tuning, no leakage)
- **The test partition is sealed** — never seen during tuning
- **Accession-overlap filter for external validation** — any protein
  appearing in both training and an external set (PhaSePro, LLPSDB2)
  is removed; AUC is reported only on the non-overlapping subset
- Detailed overlap counts written to
  `products/<product>/leakage_report.json`

The headline numbers are in the "Key numbers" table at the top. Both
products live independently under `products/` and can be used
separately.

## Citation

This repository is the updated successor to the Tingting Li Lab's
earlier PhaSePred work. If you use this tool, please cite the original
paper:

> Chen, Z., Hou, C., Wang, L., Yu, C., Chen, T., Shen, B., Hou, Y.,
> Li, P., Li, T. (2022). Screening membraneless organelle participants
> with machine-learning models that integrate multimodal features.
> *Proceedings of the National Academy of Sciences* 119(24), e2115369119.
> https://doi.org/10.1073/pnas.2115369119

External feature tools and their citations are listed in
[docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md).

## License

PhaSePred (the code, trained models, and documentation in this repository)
is released under the **MIT License**. See [LICENSE](LICENSE) at the repo
root for the full text.

> Copyright (c) 2026 Tingting Li Lab, Department of Biochemistry and Molecular Biology,
> School of Basic Medical Sciences, Peking University.

Vendored third-party components (PLAAC, PScore, SEG) keep their upstream
licenses; see [docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md)
for the full audit. User-installed tools (ESpritz, IUPred3, DeepCoil) are
governed by their own upstream terms.
