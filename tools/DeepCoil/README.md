# DeepCoil — tool package

Self-contained DeepCoil contract (PhaSePred S3 tool standardization).

## What is DeepCoil

DeepCoil (labstructbioinf) predicts coiled-coil domains from sequence.
PhaSePred's v2022 `DeepCoil` feature is the **binary indicator** `1.0 if
max(raw_cc) >= 0.82 else 0.0` computed from the per-residue `raw_cc`
probabilities (FEATURES.md §8, matching the paper's 0.82 threshold).

- Upstream: <https://github.com/labstructbioinf/DeepCoil>, PyPI `deepcoil`
- License: **none declared** — default all-rights-reserved; **not
  redistributable**. Consumed via `pip install deepcoil==2.0.2` inside an
  isolated environment (not vendored). See `docs/TOOL_LICENSES.md` §6.
- Reference: Ludwiczak J., Winski A., Szczepaniak K., Alva V.,
  Dunin-Horkawicz S. (2019) "DeepCoil — fast and accurate prediction of
  coiled-coil domains in protein sequences" *Bioinformatics*.

## Package layout

```
tools/DeepCoil/
├── run            entrypoint: ./run [-i in.fasta -out_path OUT -n_cpu N]
├── manifest.toml  contract (kind=env-prefix, license=none, NOT redistributable)
├── install.sh     builds the isolated conda env (idempotent; --check)
└── README.md      this file
```

No entity is vendored — the tool lives in an isolated Python 3.8 conda env
at `.external_envs/deepcoil` (gitignored).

## Install

```bash
bash tools/DeepCoil/install.sh          # conda-create + pip install deepcoil==2.0.2
bash tools/DeepCoil/install.sh --check  # verify env/bin/deepcoil present
```

Overrides: `DEEPCOIL_ENV_PREFIX=/path/to/env DEEPCOIL_PYTHON_VERSION=3.8
DEEPCOIL_VERSION=2.0.2`. The env prefix declared in
`manifest.toml` `[runtime]` is `.external_envs/deepcoil`; the manifest
`env_prefix` and the wrapper's `DEEPCOIL_ENV_PREFIX` override agree.

## Usage

```bash
DEEPCOIL_ENV_PREFIX="$PWD/.external_envs/deepcoil" \
tools/DeepCoil/run -i in.fasta -out_path /tmp/dc_out -n_cpu 4
tools/DeepCoil/run --check
```

DeepCoil sanitizes output filenames (alnum + `_`). Writes
`<safe_id>.out` tab-separated files (`aa`, `cc`, `raw_cc`, `prob_a`,
`prob_d`) into `--out_path` (declared output format in `manifest.toml`:
`deepcoil`).

## Attribution

Reference the DeepCoil meta-server page / Ludwiczak et al. (2019) when
publishing scores produced by this tool.