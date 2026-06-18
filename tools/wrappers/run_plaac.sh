#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
jar="$repo_root/tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar"

if [[ ! -f "$jar" ]]; then
  cat >&2 <<EOF
PLAAC jar not found: $jar

Expected local install:
  unzip tools/per-tool/PLAAC/plaac-master.zip -d tools/per-tool/PLAAC
EOF
  exit 1
fi

exec java -jar "$jar" "$@"
