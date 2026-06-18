#!/usr/bin/env bash
# Convenience: install or verify every external tool PhaSePred needs.
#
# Tools that are vendored under permissive licenses (PLAAC, PScore,
# catGRANULE) just work after `git clone`. Tools that need a local build
# (SEG) or a user-supplied download (IUPred3, ESpritz, DeepCoil) are
# handled here.
#
# The script is idempotent: each per-tool installer no-ops when its
# target is already in place.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
install_dir="$repo_root/tools/install"

echo "==> Building SEG ..."
bash "$install_dir/install_seg.sh"

echo "==> DeepCoil conda env ..."
if command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
  bash "$install_dir/install_deepcoil_env.sh"
else
  cat <<'EOF'
SKIPPED: neither conda nor mamba was found on PATH.
DeepCoil requires an isolated Python 3.7-3.8 env (TensorFlow 2.3 era).
Either:
  - install miniconda / mambaforge, then rerun this script
  - or skip DeepCoil for now (predict will report DeepCoil:0/N missing)
EOF
fi

echo "==> ESpritz (academic-license, user-supplied archive) ..."
if [[ -f "$repo_root/tools/per-tool/ESpritz/espritz.zip" ]] && [[ ! -f "$repo_root/tools/per-tool/ESpritz/espritz/espritz.pl" ]]; then
  bash "$install_dir/install_espritz.sh"
elif [[ -f "$repo_root/tools/per-tool/ESpritz/espritz/espritz.pl" ]]; then
  echo "Already extracted at tools/per-tool/ESpritz/espritz/. Skipping."
else
  echo "SKIPPED: tools/per-tool/ESpritz/espritz.zip not present."
  echo "Place the archive there and rerun, or set PHASEPRED_ESPRITZ_DIR."
fi

echo "==> IUPred3 (academic-license, user-supplied archive) ..."
if [[ -f "$repo_root/tools/per-tool/IUPred3/iupred3.tar.gz" ]] && [[ ! -f "$repo_root/tools/per-tool/IUPred3/iupred3/iupred3_lib.py" ]]; then
  bash "$install_dir/install_iupred3.sh"
elif [[ -f "$repo_root/tools/per-tool/IUPred3/iupred3/iupred3_lib.py" ]]; then
  echo "Already extracted at tools/per-tool/IUPred3/iupred3/. Skipping."
else
  echo "SKIPPED: tools/per-tool/IUPred3/iupred3.tar.gz not present."
  echo "Download from https://iupred3.elte.hu/ and rerun, or set PHASEPRED_IUPRED3_DIR."
fi

echo
echo "==> Verifying all tools resolve correctly ..."
cd "$repo_root"
uv run phasepred check-tools
