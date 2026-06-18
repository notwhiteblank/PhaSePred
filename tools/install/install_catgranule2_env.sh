#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tool_dir="$repo_root/tools/per-tool/catGRANULE_v2/catGRANULE2.0-1.0.0"
env_prefix="${CATGRANULE2_ENV_PREFIX:-$repo_root/.external_envs/catgranule2}"
python_version="${CATGRANULE2_PYTHON_VERSION:-3.10}"

if [[ ! -f "$tool_dir/requirements.txt" ]]; then
  cat >&2 <<EOF
catGRANULE2 requirements not found: $tool_dir/requirements.txt

Expected local archive extraction:
  tar -xzf tools/per-tool/catGRANULE_v2/catGRANULE2.0-1.0.0.tar.gz -C tools/per-tool/catGRANULE_v2
EOF
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  cat >&2 <<'EOF'
conda is required to create the isolated catGRANULE2 environment.
Install Miniforge/Miniconda or put conda on PATH, then rerun this script.
EOF
  exit 1
fi

mkdir -p "$(dirname "$env_prefix")"

if [[ ! -x "$env_prefix/bin/python" ]]; then
  conda create --yes --prefix "$env_prefix" "python=$python_version" pip
fi

"$env_prefix/bin/python" -m pip install --upgrade "pip<25" "setuptools<81" "wheel<0.46"
SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True \
  "$env_prefix/bin/python" -m pip install -r "$tool_dir/requirements.txt"

# training_catGRANULE2.py imports catboost, but upstream requirements.txt omits it.
"$env_prefix/bin/python" -m pip install "catboost>=1.2,<1.3"
"$env_prefix/bin/python" -m pip install "setuptools<81" "wheel<0.46"

"$env_prefix/bin/python" - <<'PY'
import Bio
import joblib
import numpy
import pandas
import scipy
import sklearn
import xgboost

print("catGRANULE2 environment ok")
PY
