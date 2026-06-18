# DeepCoil

DeepCoil upstream:

- GitHub: `https://github.com/labstructbioinf/DeepCoil`
- PyPI/install command advertised upstream: `pip install deepcoil`

## Local Decision

Do not install DeepCoil into the project `uv` environment.

As checked on 2026-05-07, upstream declares `python>=3.7,<3.9` and pins legacy
machine-learning dependencies, including `allennlp==0.9.0`, TensorFlow 2.3-era
requirements, `pandas==1.3.0`, and `biopython==1.79`. This conflicts with this
repository's Python 3.12 project environment and current scientific stack.

Treat DeepCoil as an external feature tool, not a project Python dependency.

## Isolated Environment

Create the local conda prefix environment with:

```bash
tools/install/install_deepcoil_env.sh
```

Default environment prefix:

```text
.external_envs/deepcoil
```

Override the prefix, Python version, or package version if needed:

```bash
DEEPCOIL_ENV_PREFIX=/path/to/env DEEPCOIL_PYTHON_VERSION=3.8 DEEPCOIL_VERSION=2.0.2 \
  tools/install/install_deepcoil_env.sh
```

Run the CLI through the wrapper:

```bash
tools/wrappers/run_deepcoil.sh -h
```

## Intended Integration

For the reconstruction pipeline, run DeepCoil from the isolated Python 3.8
environment or a dedicated CPU-only container. The main PhaSePred code invokes
it as a command-line tool on FASTA input and parses the generated output files,
rather than importing `deepcoil`.

Use the batch runner for recomputed PhaSePred sequence tables:

```bash
DEEPCOIL_BATCH_SIZE=25 DEEPCOIL_N_CPU=4 the project DeepCoil batch runner
```

Generated DeepCoil feature outputs belong under `data/interim/` and
`runs/deepcoil_scores/`, not under this tool archive.

## Local Status on 2026-05-11

- CLI help works through `tools/wrappers/run_deepcoil.sh`.
- Smoke scoring works; a 50-sequence CPU run took about 2 minutes 22 seconds
  with `DEEPCOIL_N_CPU=4`.
- `The project DeepCoil aggregation script` aggregates per-residue output as:
  `DeepCoil` = mean `raw_cc`, `DeepCoil_max` = max `raw_cc`, and
  `DeepCoil_fraction_cc` = fraction of residues with sharpened `cc > 0`.
- Full scoring over 60,711 recomputed sequence rows is expected to take roughly
  48 hours at the smoke-test rate and was not completed during setup.

## Open Work

- Decide whether to mirror the upstream source/archive locally.
- Run or resume full 60,711-sequence scoring when compute time is available.
- Decide fallback policy for `legacy` versus `modern` coiled-coil features.
