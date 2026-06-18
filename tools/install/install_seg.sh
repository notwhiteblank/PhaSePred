#!/usr/bin/env bash
# Build the vendored NCBI SEG source into a usable binary.
#
# Vendored at tools/per-tool/SEG/ (NCBI public-domain convention; see
# docs/THIRD_PARTY_LICENSES.md for the audit). Already shipped with the
# repository — this script just runs `make`.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
seg_dir="$repo_root/tools/per-tool/SEG"

if [[ ! -d "$seg_dir" ]]; then
  echo "SEG source dir missing: $seg_dir" >&2
  exit 1
fi

if [[ ! -f "$seg_dir/seg.c" ]]; then
  echo "SEG source files missing under $seg_dir. Did the vendored tree get deleted?" >&2
  exit 1
fi

if [[ -x "$seg_dir/seg" ]]; then
  echo "SEG already built at $seg_dir/seg — skipping. (Force a rebuild with:"
  echo "  cd $seg_dir && make clean && make )"
  exit 0
fi

echo "Building SEG in $seg_dir ..."
(cd "$seg_dir" && make)

if [[ ! -x "$seg_dir/seg" ]]; then
  echo "Build did not produce $seg_dir/seg" >&2
  exit 1
fi

echo "Done. Verify with:"
echo "  uv run phasepred check-tools"
