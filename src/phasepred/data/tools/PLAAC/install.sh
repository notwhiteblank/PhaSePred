#!/usr/bin/env bash
# PLAAC installer — the prebuilt jar is vendored; this verifies/extracts it.
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
jar="$pkg_dir/plaac-master/web/bin/plaac.jar"
archive="$pkg_dir/plaac-master.zip"

if [[ "${1:-}" == "--check" ]]; then
  if [[ ! -f "$jar" ]]; then
    echo "PLAAC: jar missing ($jar). Run: bash $pkg_dir/install.sh" >&2
    exit 1
  fi
  if ! command -v java >/dev/null 2>&1; then
    echo "PLAAC: java runtime not on PATH (needs Java >= 11)" >&2
    exit 1
  fi
  echo "PLAAC: ok ($jar)"
  exit 0
fi

if [[ -f "$jar" ]]; then
  echo "PLAAC jar already in place — skipping."
else
  if [[ ! -f "$archive" ]]; then
    cat >&2 <<EOF
plaac-master.zip not found at $archive

The upstream MIT-licensed source + prebuilt jar are normally vendored inside
the repository. Restore the tree or set PHASEPRED_PLAAC_DIR to a
tools/PLAAC/ clone.
EOF
    exit 1
  fi
  echo "Extracting $archive ..."
  unzip -q -o "$archive" -d "$pkg_dir"
fi

if [[ ! -f "$jar" ]]; then
  echo "Unexpected archive layout — $jar not found after extraction." >&2
  exit 1
fi

echo "Done. Verify with: bash $pkg_dir/install.sh --check"