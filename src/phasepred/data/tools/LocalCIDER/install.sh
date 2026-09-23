#!/usr/bin/env bash
# LocalCIDER installer — pip-installs localcider from PyPI (GPL-2.0).
# Consumed as a dependency, never vendored into the repository.
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Mirrors LocalCIDER/run's resolver (Ruling 1): $PHASEPRED_PYTHON first, then
# common interpreter names on PATH. Keep the two lists in step so --check
# cannot disagree with the runner. This wheel-bundled copy (D40) deliberately
# omits the checkout's repo-relative virtualenv candidates, whose inferred
# repository root would be the package's own data directory here rather than
# the checkout.
python_candidates() {
  local cand
  if [[ -n "${PHASEPRED_PYTHON:-}" ]]; then
    echo "$PHASEPRED_PYTHON"
  fi
  for cand in python3.13 python3.12 python3 python; do
    if command -v "$cand" >/dev/null 2>&1; then
      echo "$cand"
    fi
  done
}

resolve_python() {
  local cand
  while IFS= read -r cand; do
    [[ -z "$cand" ]] && continue
    if command -v "$cand" >/dev/null 2>&1 || [[ -x "$cand" ]]; then
      echo "$cand"
      return 0
    fi
  done < <(python_candidates)
  return 1
}

if [[ "${1:-}" == "--check" ]]; then
  while IFS= read -r cand; do
    [[ -z "$cand" ]] && continue
    if "$cand" -c "import localcider" >/dev/null 2>&1; then
      echo "LocalCIDER: ok via $cand"
      exit 0
    fi
  done < <(python_candidates)
  echo "LocalCIDER: not importable. pip install localcider==0.1.21 (GPL-2.0, PyPI)" >&2
  exit 1
fi

py="$(resolve_python)" || {
  echo "No python3/python interpreter found on PATH." >&2
  exit 1
}

if "$py" -c "import localcider" >/dev/null 2>&1; then
  echo "localcider already importable via $py — skipping."
else
  echo "Installing localcider==0.1.21 into $py ..."
  "$py" -m pip install "localcider==0.1.21"
fi

echo "Done. Verify with: bash $pkg_dir/install.sh --check"