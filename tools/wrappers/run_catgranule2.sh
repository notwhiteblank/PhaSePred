#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tool_dir="$repo_root/tools/per-tool/catGRANULE_v2/catGRANULE2.0-1.0.0"
env_prefix="${CATGRANULE2_ENV_PREFIX:-$repo_root/.external_envs/catgranule2}"
entrypoint="$tool_dir/compute_profiles_and_predictions.py"

if [[ ! -x "$env_prefix/bin/python" ]]; then
  cat >&2 <<EOF
catGRANULE2 environment not found: $env_prefix

Create it with:
  tools/install/install_catgranule2_env.sh
EOF
  exit 1
fi

if [[ ! -f "$entrypoint" ]]; then
  cat >&2 <<EOF
catGRANULE2 entrypoint not found: $entrypoint

Expected local archive extraction:
  tar -xzf tools/per-tool/catGRANULE_v2/catGRANULE2.0-1.0.0.tar.gz -C tools/per-tool/catGRANULE_v2
EOF
  exit 1
fi

if [[ ! -e "$tool_dir/src/ChemicalPhysicalScales_Py_dictionary" ]]; then
  ln -s ../ChemicalPhysicalScales_Py_dictionary "$tool_dir/src/ChemicalPhysicalScales_Py_dictionary"
fi

if [[ "${1:-}" == "--smoke-test" ]]; then
  cd "$tool_dir"
  exec "$env_prefix/bin/python" - <<'PY'
import os
import joblib
import numpy
import pandas
import sklearn
import xgboost

from compute_profiles_and_predictions import correct_order_columns

model_path = "src/TRAINED_MODELS/ONLY_PHYSCHEM/RandomForest/gridsearchCV_Object.pkl"
scale_dir = "src/ChemicalPhysicalScales_Py_dictionary"
missing = [
    name
    for name in correct_order_columns[:82]
    if name not in {"charge", "fg", "rg"} and not os.path.exists(f"{scale_dir}/{name}.json")
]

joblib.load(model_path)
print("catGRANULE2 imports and bundled RandomForest model ok")
if missing:
    print(f"warning: {len(missing)} required physchem scale JSON files are missing from this release")
PY
fi

cd "$tool_dir"
exec "$env_prefix/bin/python" "$@"
