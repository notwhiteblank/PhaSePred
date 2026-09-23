#!/usr/bin/env bash
# SEG installer — builds the vendored NCBI SEG binary (idempotent).
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${1:-}" == "--check" ]]; then
  if [[ -x "$pkg_dir/seg" ]]; then
    echo "SEG: ok ($pkg_dir/seg)"
    exit 0
  fi
  echo "SEG: not built ($pkg_dir/seg). Run: bash $pkg_dir/install.sh" >&2
  exit 1
fi

if [[ ! -d "$pkg_dir" ]]; then
  echo "SEG source dir missing: $pkg_dir" >&2
  exit 1
fi

if [[ ! -f "$pkg_dir/seg.c" ]]; then
  echo "SEG source files missing under $pkg_dir. Did the vendored tree get deleted?" >&2
  exit 1
fi

if [[ -x "$pkg_dir/seg" ]]; then
  echo "SEG already built at $pkg_dir/seg — skipping. (Force a rebuild with:"
  echo "  cd $pkg_dir && make clean && make )"
  exit 0
fi

echo "Building SEG in $pkg_dir ..."
(cd "$pkg_dir" && make)

if [[ ! -x "$pkg_dir/seg" ]]; then
  echo "Build did not produce $pkg_dir/seg" >&2
  exit 1
fi

echo "Done. Verify with:"
echo "  bash $pkg_dir/install.sh --check"