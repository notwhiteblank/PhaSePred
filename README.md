# PhaSePred

> Predict the phase-separation propensity of a protein from its amino-acid sequence.
> Open implementation of [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119).
> Maintained by the Tingting Li Lab, Peking University.

[中文](README-zh.md) · [![ci](https://github.com/NotWhiteBlank/PhaSePred/actions/workflows/ci.yml/badge.svg)](https://github.com/NotWhiteBlank/PhaSePred/actions/workflows/ci.yml)

## What PhaSePred is

**Phase separation** is the process by which proteins spontaneously condense into droplet-like compartments inside a cell — the way oil forms droplets in water. These "membraneless organelles" (stress granules, the nucleolus, Cajal bodies) organise biochemistry without a surrounding membrane. Dysregulated phase separation is closely linked to neurodegenerative disease (ALS, Alzheimer's) and to cancer.

PhaSePred predicts a protein's propensity to undergo phase separation from its sequence, and distinguishes two mechanisms across four modes:

| Mode | Meaning | Scope |
|------|---------|-------|
| **SaPS** | Self-driven phase separation — the protein condenses on its own | any species (8 features) |
| **PdPS** | Partner-dependent phase separation — requires a binding partner | any species (8 features) |
| **hSaPS** | Human self-driven phase separation | human proteins (10 features) |
| **hPdPS** | Human partner-dependent phase separation | human proteins (10 features) |

The output is a score between 0 and 1; higher means more likely to belong to that mode.

The ten features, of which the 8-feature modes use all but the last two:

| Feature | Meaning | Computed by | 8f | 10f |
|---------|---------|-------------|:--:|:---:|
| **Hydropathy** | per-residue normalised Kyte-Doolittle hydropathy mean | LocalCIDER | ✓ | ✓ |
| **FCR** | fraction of charged residues | LocalCIDER | ✓ | ✓ |
| **IDR** | intrinsically disordered region fraction | ESpritz (DisProt model `D`, `sw 0`, 5% FPR) | ✓ | ✓ |
| **LCR** | low-complexity region fraction | SEG | ✓ | ✓ |
| **PScore** | amino-acid composition complexity score | PScore | ✓ | ✓ |
| **PLAAC** | prion-like domain NLLR score | PLAAC | ✓ | ✓ |
| **catGRANULE** | granule score after Bolognesi et al. 2016 | the `catgranule` package | ✓ | ✓ |
| **DeepCoil** | coiled-coil, binarised at the paper's 0.82 threshold | DeepCoil | ✓ | ✓ |
| **Phos freq** | phosphorylation site frequency | PhosphoSitePlus | | ✓ |
| **DeepPhase** | DeepPhase phase-separation score | DeepPhase (ships with the package) | | ✓ |

## Quick start

### Install

The package needs Python `>=3.12`. Pick one of four routes:

```bash
pixi install                                                        # pixi
uv venv --python 3.12 && uv pip install .                           # uv
conda env create -f environment.yml && conda activate phasepred     # conda
python -m venv .venv && .venv/bin/pip install -e . -e packages/catgranule   # plain venv
```

Install the external tools separately:

```bash
bash tools/SEG/install.sh
bash tools/PLAAC/install.sh
bash tools/PScore/install.sh
bash tools/ESpritz/install.sh
bash tools/DeepCoil/install.sh
bash tools/LocalCIDER/install.sh
bash tools/PhosphoSitePlus/install.sh   # only needed for hSaPS / hPdPS
```

Then check:

```bash
phasepred check-tools          # table, with an install hint on every MISSING row
phasepred check-tools --strict # exits 1 when a required tool is missing; use this in scripts and CI
```

`SEG`, `PLAAC`, `LocalCIDER`, `catGRANULE` and `DeepPhase` ship with the package and work immediately. `PScore`, `ESpritz`, `DeepCoil` and `PhosphoSitePlus` need the installers above.

### Run

```bash
# 8-feature modes (any species)
phasepred predict --fasta proteins.fasta --mode SaPS --output scores.csv
phasepred predict --fasta proteins.fasta --mode PdPS --output scores.csv

# 10-feature modes (human proteins)
phasepred predict --fasta proteins.fasta --mode hSaPS --output scores.csv
phasepred predict --fasta proteins.fasta --mode hPdPS --output scores.csv

# UniProt ID input
phasepred predict --ids "P35637,Q9Y2W1" --mode SaPS --output scores.csv    # inline
phasepred predict --ids my_ids.txt --mode SaPS --output scores.csv         # file, one per line
cat my_ids.txt | phasepred predict --ids - --mode SaPS --output scores.csv # stdin

# FASTA and IDs can be combined
phasepred predict --fasta my.fasta --ids "P35637" --mode SaPS --output scores.csv

# features only, no prediction
phasepred features-from-fasta --input my.fasta --output features.csv
```

The output CSV carries the numeric value of every feature alongside `score`.

## Reproduction

The four training tables are in `data/processed/`: a serialisation of the paper's Dataset S2/S3 feature columns (base/human × train/test).

Training from scratch and AUC validation are one command each, from the repository root:

```bash
./train.sh      # retrain the 40 models and byte-compare against the shipped artifacts
./validate.sh   # print the four-mode AUC table against the paper
```

`train.sh` writes the 40 models to `runs/retrain/`, sha256-compares each against `src/phasepred/data/models/`, and ends with an `N/40 byte-identical` line. Exit codes: `0` all identical, `1` an artifact differs, `2` preflight failure. Byte-exact reproduction requires `xgboost>=3.2,<3.3`.

`validate.sh` runs three layers (paper AUC, leakage, artifact consistency), prints the comparison table to the screen, and writes its report to `runs/validation/`:

```
MODE        PAPER      MEASURED      DELTA  VERDICT
SaPS-8      0.862      0.860729    -0.0013  pass
PdPS-8      0.739      0.737884    -0.0011  pass
hSaPS-10    0.924      0.917515    -0.0065  pass
hPdPS-10    0.827      0.816171    -0.0108  fail
```

Exit codes: `0` no hard gate violated, `1` a hard gate failed, `2` input or usage error. Metrics and splits land in `products/A_paper_split_recomputed/`.

## Citation

> Chen, Z., Hou, C., Wang, L., Yu, C., Chen, T., Shen, B., Hou, Y.,
> Li, P., Li, T. (2022). Screening membraneless organelle participants
> with machine-learning models that integrate multimodal features.
> *Proceedings of the National Academy of Sciences* 119(24), e2115369119.
> https://doi.org/10.1073/pnas.2115369119

The catGRANULE scoring implementation derives from Bolognesi et al. 2016 (*Cell Reports* 16:222-231).

## License

This repository (code, models, documentation) is released under the **MIT License** — see [LICENSE](LICENSE).

> Copyright (c) 2026 Tingting Li Lab, Department of Biochemistry and Molecular Biology,
> School of Basic Medical Sciences, Peking University.

Third-party components and their licences are listed in [NOTICE](NOTICE).
