# Python Environment

This project uses [`uv`](https://docs.astral.sh/uv/) for Python environment
and lockfile management.

## Python Version

Pinned to Python 3.12 through `.python-version` and `pyproject.toml`.

## Setup

```bash
uv sync
```

If Python downloads are unreliable, use a locally available interpreter:

```bash
uv sync --no-python-downloads
```

Run Python inside the project environment:

```bash
uv run python --version
```

Quick dependency smoke check:

```bash
uv run python - <<'PY'
import Bio
import localcider
import openpyxl
import pandas
import sklearn
import xgboost

print("environment ok")
PY
```

## localCIDER

`localcider` is a Python dependency installed in the `uv` environment. It
calculates sequence properties relevant to PhaSePred features, including
hydropathy, fraction of charged residues, net charge, kappa, and
isoelectric point.

## Dependency Groups

Default dependencies cover data parsing, Excel handling, sequence
processing, modeling, and CLI development. Development dependencies live
in the `dev` group:

- `pytest`
- `ruff`
- `mypy`

## External Non-Python Tools

The `uv` environment intentionally does not install the legacy
command-line feature tools. They run through wrappers under
`tools/wrappers/` (resolved via `src/phasepred/tool_paths.py`) so the
Python environment stays small and the feature stack stays version-pinned
to whatever the user has on disk.

See [tools/README.md](../tools/README.md) for the resolution mechanics
(env var → vendored → PATH → error) and
[docs/THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for the license
audit. The short version:

| Tool | Vendored in repo? | Install | Wrapper |
|---|---|---|---|
| SEG | source | `bash tools/install/install_seg.sh` | (used directly by `tool_paths.find_seg`) |
| PScore | source | already vendored | (used directly) |
| PLAAC | prebuilt jar | already vendored; needs Java 17+ | `tools/wrappers/run_plaac.sh` |
| catGRANULE v1 | reimplemented in `src/phasepred/catgranule_v1.py` | n/a | n/a |
| ESpritz | no — academic license | `bash tools/install/install_espritz.sh` | `tools/wrappers/run_espritz.sh` |
| IUPred3 | no — academic license | `bash tools/install/install_iupred3.sh` | `tools/wrappers/run_iupred3.sh` |
| DeepCoil | no — no upstream license | `bash tools/install/install_deepcoil_env.sh` (creates `.external_envs/deepcoil`) | `tools/wrappers/run_deepcoil.sh` |

After installing, verify with:

```bash
uv run phasepred check-tools
```

## Tool Isolation Policy

Keep the main `uv` environment limited to orchestration, parsing,
modeling, and CLI code. Feature tools live in isolated environments,
mostly:

- `tools/per-tool/<TOOL>/` for vendored ones,
- the system installer for `make` / `java` / `perl`,
- `.external_envs/deepcoil/` (separate conda env) for DeepCoil.

DeepCoil cannot share the `uv` environment because its upstream pins
require Python 3.7-3.8 with TensorFlow 2.3-era dependencies.

`phasepred predict` calls `tool_paths.assert_predict_requirements()`
before reading any input, so any missing tool fails up-front with a
per-tool install hint.

## Optional human-feature data

`hSaPS` and `hPdPS` predictions use two extra columns that come from
licensed third-party tables, not from inference tools:

- **DeepPhase** scores for ~6000 reference human proteins, from the
  DeepPhase paper supplement (`tableS3.xlsx`).
- **PhosphoSitePlus** phosphorylation-site dataset
  (`Phosphorylation_site_dataset.gz`), which requires registration to
  download.

Without these, `phasepred predict --mode hSaPS|hPdPS` still runs but
emits explicit `PhaSePredMissingDataWarning` entries naming the
accessions that ended up median-imputed. See the
"Installing optional human-feature data" section in the top-level
README.

Both files belong under `data/raw/external/...` (see paths in the
README) and are kept out of git history regardless of license status.
