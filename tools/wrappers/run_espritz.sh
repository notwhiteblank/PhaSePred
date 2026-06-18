#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tool_dir="$repo_root/tools/per-tool/ESpritz/espritz"
entrypoint="$tool_dir/espritz.pl"

if [[ ! -f "$entrypoint" ]]; then
  cat >&2 <<EOF
ESpritz entrypoint not found: $entrypoint

Expected local install:
  unzip tools/per-tool/ESpritz/espritz.zip -d tools/per-tool/ESpritz
  chmod +x tools/per-tool/ESpritz/espritz/espritz.pl tools/per-tool/ESpritz/espritz/bin/disbin*
EOF
  exit 1
fi

if [[ $# -ne 3 ]]; then
  cat >&2 <<EOF
Usage: tools/wrappers/run_espritz.sh WORKDIR MODEL SW

Arguments:
  WORKDIR  Directory containing one or more .fasta files.
  MODEL    X, D, N for sequence-only models; pX, pD, pN require PSI-BLAST setup.
  SW       1 for best-Sw threshold, 0 for 5% false-positive-rate threshold.

Example:
  tools/wrappers/run_espritz.sh tools/per-tool/ESpritz/espritz/example_fastas D 0
EOF
  exit 2
fi

workdir="$1"
model="$2"
sw="$3"

case "$model" in
  X|D|N) ;;
  pX|pD|pN)
    cat >&2 <<EOF
ESpritz PSI-BLAST models require additional BLAST binaries, databases, and
tools/per-tool/ESpritz/espritz/align/getAlignments.pl path configuration.
Use X, D, or N until that external database layer is prepared.
EOF
    exit 2
    ;;
  *)
    echo "Unknown ESpritz model: $model" >&2
    exit 2
    ;;
esac

case "$sw" in
  0|1) ;;
  *)
    echo "SW must be 0 or 1, got: $sw" >&2
    exit 2
    ;;
esac

if [[ ! -d "$workdir" ]]; then
  echo "ESpritz input directory not found: $workdir" >&2
  exit 1
fi

workdir_abs="$(cd "$workdir" && pwd)"

(
  cd "$tool_dir"
  exec perl "$entrypoint" "$workdir_abs" "$model" "$sw"
)
