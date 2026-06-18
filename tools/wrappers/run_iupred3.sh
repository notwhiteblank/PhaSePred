#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tool_dir="$repo_root/tools/per-tool/IUPred3/iupred3"
entrypoint="$tool_dir/iupred3.py"

if [[ ! -f "$entrypoint" ]]; then
  cat >&2 <<EOF
IUPred3 entrypoint not found: $entrypoint

Expected local install:
  tar -xzf tools/per-tool/IUPred3/iupred3.tar.gz -C tools/per-tool/IUPred3
EOF
  exit 1
fi

exec uv run python "$entrypoint" "$@"
