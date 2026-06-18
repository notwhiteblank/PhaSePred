#!/usr/bin/env bash
# IUPred3 installer — guides the user to download the academic-license
# archive and extract it to the expected vendored location.
#
# IUPred3 is distributed under an ELTE academic license that forbids
# redistribution (see docs/THIRD_PARTY_LICENSES.md), so this script only
# moves a tarball you have already downloaded. If you'd rather install
# IUPred3 elsewhere, set PHASEPRED_IUPRED3_DIR and skip this script.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
dest_dir="$repo_root/tools/per-tool/IUPred3"
archive="$dest_dir/iupred3.tar.gz"

cat <<EOF
IUPred3 cannot be redistributed under its academic license.

1. Open https://iupred3.elte.hu/ in a browser.
2. Accept the academic-use terms and download iupred3.tar.gz.
3. Place the downloaded archive at:
     $archive

Press Enter once that file is in place, or Ctrl+C to abort.
EOF
read -r _

if [[ ! -f "$archive" ]]; then
  echo "Archive not found at $archive. Aborting." >&2
  exit 1
fi

echo "Extracting $archive into $dest_dir ..."
tar -xzf "$archive" -C "$dest_dir"

if [[ ! -d "$dest_dir/iupred3" ]] || [[ ! -f "$dest_dir/iupred3/iupred3_lib.py" ]]; then
  echo "Unexpected archive layout — iupred3/iupred3_lib.py not found." >&2
  exit 1
fi

echo "Done. Verify with:"
echo "  uv run phasepred check-tools"
