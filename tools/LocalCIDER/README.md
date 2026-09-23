# LocalCIDER — tool package

Self-contained LocalCIDER contract (PhaSePred S3 tool standardization).

## What is LocalCIDER

LocalCIDER (Pappu lab) computes sequence-composition/patterning properties
of intrinsically disordered proteins. PhaSePred uses its Uversky hydropathy
and FCR definitions for the native `Hydropathy`/`FCR` 8-feature columns
(in-process import in `src/phasepred/features.py`). This contract exposes
the same library as a standalone CLI (`kind=python-library`), producing a
core sequence-properties JSON:

`length`, `uversky_hydropathy`, `fcr`, `ncpr`, `kappa`, `scd`,
`mean_net_charge`, `aromatic_fraction`.

- Upstream: <https://github.com/holehouse-lab/localcider> (PyPI `localcider`)
- License: **GPL-2.0** (0.1.21 dist-info LICENSE.txt). Redistributable via
  PyPI; **consumed through pip, never vendored** into this repository
  (avoids GPL copyleft contamination). See `docs/TOOL_LICENSES.md` §7.
- Reference: Holehouse A.S., Das R.K., Ahad J.N., Richardson M.O.C.,
  Pappu R.V. (2017) CIDER: Resources to analyze sequence-ensemble
  relationships of intrinsically disordered proteins. *Biophysical
  Journal* 113:1484-1492.

## Package layout

```
tools/LocalCIDER/
├── run            entrypoint: ./run <in.fasta>  → JSON on stdout
├── manifest.toml  contract (kind=python-library, license=GPL-2.0)
├── install.sh     pip install localcider==0.1.21 (idempotent; --check)
└── README.md      this file
```

## Install

```bash
bash tools/LocalCIDER/install.sh          # pip install localcider==0.1.21
bash tools/LocalCIDER/install.sh --check  # verify importable
```

**Resolving the interpreter is the caller's job.** A bash script cannot know
its caller's `sys.executable`, so `PHASEPRED_PYTHON` is the recommended and
authoritative choice: when set, it is tried first and used if it can import
`localcider`. Every subprocess call in `src/phasepred/` (13 sites) passes
`PHASEPRED_PYTHON=sys.executable` through `contract_env()`, so the contract
runs under the same interpreter as `predict`. Without the variable,
`resolve_python()` falls back to a widened list — `python3.13` → `python3.12`
→ `python3` → `python`, then the repo-relative `.venv/bin/python` and
`.pixi/envs/default/bin/python` — requiring an importable `localcider` at each
candidate. `install.sh --check` mirrors this list.

## Usage

```bash
tools/LocalCIDER/run tests/fixtures/sequences/Q08211.fasta
tools/LocalCIDER/run --check
```

Standard output is a JSON array (declared output format in
`manifest.toml`: `json`). On a missing `localcider`, the wrapper prints the
install hint and exits non-zero.

`src/phasepred/features.py` keeps importing `localcider` in-process; this
`run` entry is for standalone use and discoverability.