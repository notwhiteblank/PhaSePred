# catGRANULE v2

Local archive:

- `catGRANULE2.0-1.0.0.tar.gz`
- extracted tree: `catGRANULE2.0-1.0.0/`

This is the practical local substitute for paper-era `catGRANULE` while the
older source remains unavailable. Treat scores from this release as a modern
tool feature, not as the original PhaSePred paper `catGRANULE` feature.

## Isolated Environment

Do not install this dependency set into the main project `uv` environment.
Create the local conda prefix environment with:

```bash
tools/install/install_catgranule2_env.sh
```

Default environment prefix:

```text
.external_envs/catgranule2
```

Override the prefix or Python version if needed:

```bash
CATGRANULE2_ENV_PREFIX=/path/to/env CATGRANULE2_PYTHON_VERSION=3.10 \
  tools/install/install_catgranule2_env.sh
```

## Wrapper

Run Python scripts inside the extracted upstream directory with:

```bash
tools/wrappers/run_catgranule2.sh SCRIPT_OR_ARGS...
```

Smoke test the bundled TDP43 example:

```bash
tools/wrappers/run_catgranule2.sh --smoke-test
```

The wrapper changes into the upstream source directory before execution because
the released code uses relative paths such as `./src/TRAINED_MODELS/`.
It also creates a compatibility symlink from
`src/ChemicalPhysicalScales_Py_dictionary` to the top-level scale dictionary if
needed.

## Notes

- The upstream `requirements.txt` omits `catboost`, but
  `training_catGRANULE2.py` imports it. The installer adds `catboost>=1.2,<1.3`
  for training compatibility.
- The local `catGRANULE2.0-1.0.0` release and the upstream GitHub `v1.0.0`
  tree both contain 53 chemical/physical scale JSON files, while
  `compute_profiles_and_predictions.py` references 82 scale-derived feature
  names. Environment and model import smoke tests pass, but direct FASTA
  prediction currently fails until those missing scale definitions are
  recovered or the released feature set is reconciled.
- The release includes `src/catG2_scores_human_proteome.csv.zip`. Use
  `the project catGRANULE2 export utility` to export those
  precomputed human-proteome scores as the current practical `catGRANULE`
  feature table.
- Generated feature tables belong under a future derived/interim data
  lifecycle directory, not inside `tools/per-tool/catGRANULE_v2/`.
