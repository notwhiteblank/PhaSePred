#!/usr/bin/env bash
# ESpritz installer — guides the user to drop the academic-license archive
# in place and extracts it to the expected vendored location.
#
# ESpritz is distributed under the Tosatto lab academic license that
# forbids redistribution (see docs/THIRD_PARTY_LICENSES.md). If you'd
# rather install ESpritz elsewhere, set PHASEPRED_ESPRITZ_DIR and skip
# this script.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
dest_dir="$repo_root/tools/per-tool/ESpritz"
archive="$dest_dir/espritz.zip"

cat <<EOF
ESpritz cannot be redistributed under its academic license.

1. Request access from the Tosatto lab (contact details in the LICENSE
   file shipped inside the archive). Their public website is at
   https://protein.bio.unipd.it/espritz/.
2. Place the downloaded espritz.zip at:
     $archive
3. Make sure perl is installed and on PATH.

Press Enter once the archive is in place, or Ctrl+C to abort.
EOF
read -r _

if [[ ! -f "$archive" ]]; then
  echo "Archive not found at $archive. Aborting." >&2
  exit 1
fi

echo "Extracting $archive into $dest_dir ..."
unzip -q -o "$archive" -d "$dest_dir"

if [[ ! -d "$dest_dir/espritz" ]] || [[ ! -f "$dest_dir/espritz/espritz.pl" ]]; then
  echo "Unexpected archive layout — espritz/espritz.pl not found." >&2
  exit 1
fi

chmod +x "$dest_dir/espritz/espritz.pl"
chmod +x "$dest_dir/espritz/bin/"disbin* 2>/dev/null || true
chmod +x "$dest_dir/espritz/bin/blastredo" 2>/dev/null || true

echo "Done. Verify with:"
echo "  uv run phasepred check-tools"
