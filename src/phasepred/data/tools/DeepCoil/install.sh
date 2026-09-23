#!/usr/bin/env bash
# DeepCoil installer — builds an isolated Python 3.8 conda env holding
# deepcoil==2.0.2 (upstream has no LICENSE; not redistributable).
#
#   bash install.sh             build/reuse the env (fully idempotent)
#   bash install.sh --check     verify only; exit 0 installed / 1 missing
#   bash install.sh --offline   never create/download; exit 1 if unsatisfied
#
# Prefix resolution (Ruling 10):
#   $DEEPCOIL_ENV_PREFIX
#   > <user_data_root>/envs/deepcoil            (created on a fresh machine)
#   > <repo_root>/.external_envs/deepcoil       (legacy fallback, still probed)
#
# The legacy tier is load-bearing: existing checkouts already hold the 6.7 GB
# env there, and dropping it would make the DeepCoil row MISSING and break
# predict. Nothing is written into the checkout by this script.
#
# An 8 GiB free-space preflight guards a prefix that would actually be created
# (a failed conda create leaves gigabytes of residue); pass DEEPCOIL_ENV_PREFIX
# to place it on a larger filesystem.
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$pkg_dir/../.." && pwd)"
python_version="${DEEPCOIL_PYTHON_VERSION:-3.8}"
deepcoil_version="${DEEPCOIL_VERSION:-2.0.2}"

mode="install"
offline=0
for arg in "$@"; do
  case "$arg" in
    --check) mode="check" ;;
    --offline) offline=1 ;;
    -h|--help)
      sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Usage: bash $0 [--check] [--offline]" >&2
      exit 2
      ;;
  esac
done

user_root() {
  if [[ -n "${PHASEPRED_DATA_ROOT:-}" ]]; then
    printf '%s\n' "$PHASEPRED_DATA_ROOT"
  elif [[ -n "${XDG_DATA_HOME:-}" ]]; then
    printf '%s\n' "$XDG_DATA_HOME/phasepred"
  elif [[ -n "${HOME:-}" ]]; then
    printf '%s\n' "$HOME/.local/share/phasepred"
  else
    return 1
  fi
}

resolve_prefix() {
  if [[ -n "${DEEPCOIL_ENV_PREFIX:-}" ]]; then
    printf '%s\n' "$DEEPCOIL_ENV_PREFIX"
    return 0
  fi
  local user legacy
  user="$(user_root)/envs/deepcoil" || return 1
  if [[ -x "$user/bin/deepcoil" || -x "$user/bin/python" ]]; then
    printf '%s\n' "$user"
    return 0
  fi
  legacy="$repo_root/.external_envs/deepcoil"
  if [[ -x "$legacy/bin/deepcoil" || -x "$legacy/bin/python" ]]; then
    printf '%s\n' "$legacy"
    return 0
  fi
  printf '%s\n' "$user"
}

nearest_existing() {
  local p="$1"
  while [[ ! -d "$p" && "$p" != "/" ]]; do
    p="$(dirname "$p")"
  done
  printf '%s\n' "$p"
}

check_shebang() {
  local prefix="$1" exe="$prefix/bin/deepcoil" first interp
  if [[ ! -x "$exe" ]]; then
    echo "DeepCoil: $exe is missing or not executable" >&2
    return 1
  fi
  first="$(head -n1 "$exe" 2>/dev/null || true)"
  if [[ "$first" != '#!'* ]]; then
    echo "DeepCoil: $exe has no shebang; cannot verify relocatability." >&2
    return 1
  fi
  interp="${first#\#\!}"
  interp="${interp%% *}"
  case "$interp" in
    "$prefix/"*) ;;
    *)
      echo "DeepCoil: shebang points outside the prefix: $interp" >&2
      echo "  Recreate the env at this prefix rather than moving it." >&2
      return 1
      ;;
  esac
  if [[ ! -x "$interp" ]]; then
    echo "DeepCoil: shebang interpreter missing or not executable: $interp" >&2
    return 1
  fi
  return 0
}

self_check() {
  local prefix="$1"
  check_shebang "$prefix"
  if ! "$prefix/bin/deepcoil" --help >/dev/null 2>&1; then
    echo "DeepCoil: 'deepcoil --help' failed; the env is incomplete." >&2
    return 1
  fi
}

if ! prefix="$(resolve_prefix)"; then
  echo "DeepCoil: no data root (set HOME, PHASEPRED_DATA_ROOT, or DEEPCOIL_ENV_PREFIX)." >&2
  exit 1
fi

if [[ "$mode" == "check" ]]; then
  if [[ -x "$prefix/bin/deepcoil" ]]; then
    self_check "$prefix" || exit 1
    echo "DeepCoil: ok ($prefix/bin/deepcoil)"
    exit 0
  fi
  echo "DeepCoil: MISSING ($prefix/bin/deepcoil)" >&2
  echo "  Build it with: bash $pkg_dir/install.sh" >&2
  exit 1
fi

# Idempotency: a satisfied prefix is never rebuilt or re-pinned.
if [[ -x "$prefix/bin/deepcoil" ]]; then
  self_check "$prefix"
  echo "DeepCoil: already installed ($prefix) — skipping."
  exit 0
fi

if [[ "$offline" -eq 1 ]]; then
  echo "DeepCoil: --offline and the env is not installed; cannot satisfy." >&2
  echo "  Creating the env needs conda + PyPI network access." >&2
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  cat >&2 <<'EOF'
conda is required to create the isolated DeepCoil environment.
Install Miniforge/Miniconda or put conda on PATH, then rerun this script.
EOF
  exit 1
fi

# Space preflight only for a prefix that would be created. 8 GiB = 8388608 KiB.
will_create=1
[[ -x "$prefix/bin/python" ]] && will_create=0
if [[ "$will_create" -eq 1 ]]; then
  probe="$(nearest_existing "$(dirname "$prefix")")"
  avail_kb="$(df -Pk "$probe" 2>/dev/null | awk 'NR==2 {print $4}' || true)"
  if [[ -n "$avail_kb" && "$avail_kb" -lt 8388608 ]]; then
    echo "DeepCoil: not enough free space to create the env at $prefix" >&2
    echo "  $probe has $((avail_kb / 1024)) MiB free; 8192 MiB required." >&2
    echo "  Set DEEPCOIL_ENV_PREFIX to a larger filesystem and retry." >&2
    exit 1
  fi
fi

if ! mkdir -p "$(dirname "$prefix")" 2>/dev/null; then
  echo "DeepCoil: cannot create parent of $prefix (read-only?)." >&2
  echo "  Set DEEPCOIL_ENV_PREFIX to a writable path and retry." >&2
  exit 1
fi

if [[ ! -x "$prefix/bin/python" ]]; then
  conda create --yes --prefix "$prefix" "python=$python_version" pip
fi

"$prefix/bin/python" -m pip install --upgrade "pip<25" "setuptools==57.5.0" "wheel==0.32.3"

# allennlp 0.9.0 imports spacy at module import time, and spacy 2.1.9 needs
# this older build stack on modern pip/setuptools.
"$prefix/bin/python" -m pip install "Cython<3" "numpy==1.23.5"
"$prefix/bin/python" -m pip install "spacy==2.1.9" --no-build-isolation

"$prefix/bin/python" -m pip install \
  "biopython==1.79" \
  "overrides==3.1.0" \
  "pandas==1.3.0" \
  "seaborn==0.12.2" \
  "scipy==1.10.1" \
  "matplotlib==3.7.5" \
  "h5py==3.11.0" \
  "torch==1.13.1" \
  "tensorflow==2.10.1"

"$prefix/bin/python" -m pip install \
  "allennlp==0.9.0" \
  "deepcoil==$deepcoil_version" \
  --no-deps

# Imported by allennlp's package initializers even though DeepCoil only uses
# allennlp.commands.elmo.ElmoEmbedder.
"$prefix/bin/python" -m pip install \
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

self_check "$prefix"
echo "Done. Verify with: bash $pkg_dir/install.sh --check"
