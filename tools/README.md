# External Tools

PhaSePred's prediction pipeline calls into seven external feature tools.
This directory documents how they install and how PhaSePred finds them.

## Path resolution order

For every tool the path resolver in `src/phasepred/tool_paths.py` tries,
in order:

1. **Environment variable** (e.g. `PHASEPRED_SEG_BIN`,
   `PHASEPRED_ESPRITZ_DIR`).
2. **Vendored under `tools/per-tool/<TOOL>/`** in this repository.
3. **PATH** (`shutil.which`).
4. Error with the install hint string.

Run `uv run phasepred check-tools` after installing to verify every tool
is found.

## Directory layout

```
tools/
├── README.md                    this file
├── THIRD_PARTY_LICENSES.md      license audit (see docs/ too)
├── install/                     install_<tool>.sh — automated installers
├── wrappers/                    run_<tool>.sh — shell wrappers used by feature code
└── per-tool/                    per-tool README + (where licensing allows) vendored source
    ├── SEG/                     NCBI SEG (low-complexity masking)
    ├── PScore/                  eLife 2018 phase-separation predictor
    ├── PLAAC/                   prion-like amino acid composition
    ├── IUPred3/                 disorder predictor (academic license)
    ├── ESpritz/                 disorder predictor (academic license)
    ├── DeepCoil/                coiled-coil predictor (no license declared)
    └── catGRANULE_v2/           granule-formation predictor (paper-formula
                                 v1 reconstruction lives in
                                 src/phasepred/catgranule_v1.py)
```

## What's vendored vs. install-on-demand

See `THIRD_PARTY_LICENSES.md` for the full audit. Summary:

| Tool | Vendored here? | If not, install via |
|---|---|---|
| SEG | NCBI public domain (gray area) | `tools/install/install_seg.sh` (downloads + builds) |
| PScore | CC-BY 4.0, **YES** | (vendored) |
| PLAAC | MIT, **YES** | (vendored) |
| catGRANULE | paper-formula v1 reimplementation, **YES** (in `src/phasepred/catgranule_v1.py`) | n/a |
| IUPred3 | academic license forbids redistribution | `tools/install/install_iupred3.sh` (asks user to accept terms) |
| ESpritz | academic license forbids redistribution | `tools/install/install_espritz.sh` |
| DeepCoil | no license declared upstream | `tools/install/install_deepcoil_env.sh` (isolated conda env) |
| AlphaFold pLDDT | experiment-only (dropped from production) | see `experiments/2026-05-19_alphafold_plddt_ablation/` |

## After installing

```bash
uv run phasepred check-tools
```

Example output:

```
TOOL          STATUS  SOURCE                       PATH
SEG           OK      vendored                     tools/per-tool/SEG/seg
PScore        OK      vendored                     tools/per-tool/PScore/SourceCodeS2
PLAAC         OK      vendored                     tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar
IUPred3       OK      env:PHASEPRED_IUPRED3_DIR    /opt/iupred3
ESpritz       OK      env:PHASEPRED_ESPRITZ_DIR    /opt/espritz
DeepCoil      OK      vendored (conda)             .external_envs/deepcoil
catGRANULE    OK      built-in                     src/phasepred/catgranule_v1.py
```

`predict` aborts up-front if anything is MISSING.
