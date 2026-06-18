#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${DEEPCOIL_ENV_PREFIX:-$repo_root/.external_envs/deepcoil}"
python_version="${DEEPCOIL_PYTHON_VERSION:-3.8}"
deepcoil_version="${DEEPCOIL_VERSION:-2.0.2}"

if ! command -v conda >/dev/null 2>&1; then
  cat >&2 <<'EOF'
conda is required to create the isolated DeepCoil environment.
Install Miniforge/Miniconda or put conda on PATH, then rerun this script.
EOF
  exit 1
fi

mkdir -p "$(dirname "$env_prefix")"

if [[ ! -x "$env_prefix/bin/python" ]]; then
  conda create --yes --prefix "$env_prefix" "python=$python_version" pip
fi

"$env_prefix/bin/python" -m pip install --upgrade "pip<25" "setuptools==57.5.0" "wheel==0.32.3"

# allennlp 0.9.0 imports spacy at module import time, and spacy 2.1.9 needs
# this older build stack on modern pip/setuptools.
"$env_prefix/bin/python" -m pip install "Cython<3" "numpy==1.23.5"
"$env_prefix/bin/python" -m pip install "spacy==2.1.9" --no-build-isolation

"$env_prefix/bin/python" -m pip install \
  "biopython==1.79" \
  "overrides==3.1.0" \
  "pandas==1.3.0" \
  "seaborn==0.12.2" \
  "scipy==1.10.1" \
  "matplotlib==3.7.5" \
  "h5py==3.11.0" \
  "torch==1.13.1" \
  "tensorflow==2.10.1"

"$env_prefix/bin/python" -m pip install \
  "allennlp==0.9.0" \
  "deepcoil==$deepcoil_version" \
  --no-deps

# Imported by allennlp's package initializers even though DeepCoil only uses
# allennlp.commands.elmo.ElmoEmbedder.
"$env_prefix/bin/python" -m pip install \
  "boto3==1.37.38" \
  "conllu==1.3.1" \
  "editdistance==0.8.1" \
  "flaky==3.8.1" \
  "flask==3.0.3" \
  "flask-cors==5.0.0" \
  "ftfy==6.2.3" \
  "gevent==24.2.1" \
  "jsonnet==0.22.0" \
  "jsonpickle==4.1.1" \
  "nltk==3.9.1" \
  "numpydoc==1.7.0" \
  "parsimonious==0.11.0" \
  "pytest==8.3.5" \
  "pytorch-pretrained-bert==0.6.2" \
  "pytorch-transformers==1.1.0" \
  "responses==0.26.0" \
  "scikit-learn==1.3.2" \
  "sqlparse==0.5.5" \
  "tensorboardX==1.9" \
  "unidecode==1.4.0" \
  "word2number==1.1" \
  "protobuf==3.19.6"

"$env_prefix/bin/python" - <<'PY'
import deepcoil

print("DeepCoil environment ok")
PY

"$env_prefix/bin/deepcoil" -h >/dev/null
