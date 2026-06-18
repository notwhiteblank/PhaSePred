# Repository Layout

The repository is layered: a slim shared library + CLI in `src/`, two
shippable models in `products/`, and external feature tools in `tools/`.
Heavy run artifacts live locally under `runs/` and are gitignored.

## Top-level directories

```
PhaSePred/
├── src/phasepred/          shared library and CLI
├── products/               two independent product deliverables
│   ├── A_paper_split_recomputed/
│   └── B_extended_dataset/
├── tools/                  external feature tools
├── docs/                   project documentation
├── tests/                  pytest suite
├── data/                   raw + interim (gitignored)
├── runs/                   heavy run artifacts (gitignored)
└── (local-only)
    ├── PhaSePred_article&data/   paper PDFs + S2-S8 (gitignored)
    ├── support_article/          supporting paper PDFs (gitignored)
    └── .external_envs/           isolated conda envs (gitignored)
```

## `src/phasepred/` — shared library

| Module | Purpose |
|---|---|
| `cli.py` | Typer commands: `predict`, `check-tools`, `features-recomputed`, `features-from-fasta`, `features-from-split`, `features-from-paper`, `prepare-split`, `sequences-from-ids`, `train`, `predict-features` |
| `predictor.py` | FASTA → features → score pipeline; loads `products/<product>/models/` |
| `tool_paths.py` | External-tool resolver: env var → vendored → PATH → error |
| `features.py` | Native features (Hydropathy, FCR) and the schema constants `BASE_FEATURE_COLUMNS`, `HUMAN_FEATURE_COLUMNS` |
| `catgranule_v1.py` | Paper-formula catGRANULE re-implementation (built-in) |
| `iupred.py`, `espritz.py` | Feature adapters and on-disk caches |
| `legacy_features.py` | Bridges to SEG, PScore, PLAAC, PhosphoSitePlus, DeepPhase |
| `sequences.py` | UniProt FASTA fetcher + JSONL cache |
| `paper.py` | Parses paper S2/S3 supplementary workbooks (Product A) |
| `updated_data.py` | Builds the unified label table from PhaSepDB3 + PhaSePro + LLPSDB2 (Product B) |
| `models.py` | sklearn model artifact dump/load (used by the `train` CLI) |

## `products/` — shippable models

Each `<product>/` directory contains:

- `README.md` — usage instructions
- `tune.py` / `train.py` — reproducible pipeline
- `tuned_params.json` — best XGBoost params per task
- `train_accessions_<task>.tsv`, `test_accessions_<task>.tsv` — sealed
  partitions
- `leakage_report.json` — overlap counts against PhaSePro / LLPSDB2
- `metrics.json` — CV AUC, held-out test AUC, external AUCs
- `models/<task>/8f_model_{0..9}.joblib` — trained pipelines

`predict` defaults to Product B. `--product A` selects A.

## `tools/` — external tool layer

| Subdir | Contents |
|---|---|
| `README.md` | Resolver mechanics |
| `THIRD_PARTY_LICENSES.md` | License audit |
| `install/` | `install_<tool>.sh` automation scripts |
| `wrappers/` | `run_<tool>.sh` shell wrappers used by feature code |
| `per-tool/<TOOL>/` | Per-tool README; for permissively-licensed tools also the vendored source |

Tools currently vendored (per [docs/THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)):
PLAAC (MIT), PScore (CC-BY 4.0), SEG (NCBI public-domain convention).

Tools that must be installed by the user (license forbids redistribution
or upstream has no declared license): IUPred3, ESpritz, DeepCoil.

## `data/` — gitignored

```
data/
├── raw/external/          raw downloads from public databases
│   ├── phasepdb3/
│   ├── phasepro/
│   ├── llpsdb2/
│   ├── phosphositeplus/
│   ├── deepphase/
│   └── _manifests/        provenance metadata (tracked)
├── interim/               derived feature tables, label tables, caches
└── processed/             (reserved for future cleaned outputs)
```

Only `_manifests/` and the `.gitkeep` markers are tracked. Everything
else stays on disk locally.

## `runs/` — gitignored

Per-run heavy outputs (model files mid-training, plots, prediction
CSVs, log files) live locally under `runs/` and are not pushed.

## What was vendored vs left to the user

| Tool | License | In-repo? | User must do |
|---|---|---|---|
| SEG | NCBI public-domain | source vendored | `cd tools/per-tool/SEG && make` |
| PScore | CC-BY 4.0 | source vendored | nothing |
| PLAAC | MIT | prebuilt jar + cli sources vendored | install a JRE |
| catGRANULE v1 | reimplementation | built-in Python module | nothing |
| ESpritz | academic, no redistribution | README only | accept upstream terms, download zip |
| IUPred3 | academic, no redistribution | README only | accept upstream terms, download tarball |
| DeepCoil | no upstream license | README + install script | run conda env installer |
