#!/usr/bin/env bash
# ESpritz installer — unpacks the academic-license archive into the user data
# directory (idempotent).
#
#   bash install.sh             non-interactive when an archive is staged
#   bash install.sh --check     verify only; exit 0 available / 1 missing
#   bash install.sh --offline   never prompt; exit 1 if unsatisfied
#
# ESpritz is distributed under the Tosatto lab academic license (v1.1) which
# forbids redistribution (see docs/TOOL_LICENSES.md §5), so this repository
# ships no ESpritz entity and the installer never downloads. Stage
# `espritz.zip` in $PHASEPRED_VENDOR_ARCHIVE_DIR for a fully non-interactive
# install; otherwise an interactive prompt guides you through fetching it.
#
# Nothing is written into the checkout: `espritz/` and the copied `run` land
# under the user data directory ($PHASEPRED_DATA_ROOT, else
# $XDG_DATA_HOME/phasepred, else ~/.local/share/phasepred).
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
archive_name="espritz.zip"
entry_rel="espritz/espritz.pl"

mode="install"
offline=0
for arg in "$@"; do
  case "$arg" in
    --check) mode="check" ;;
    --offline) offline=1 ;;
    -h|--help)
      sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Usage: bash $0 [--check] [--offline]" >&2
      exit 2
      ;;
  esac
done

user_root() {
  if [[ -n "${PHASEPRED_DATA_ROOT:-}" ]]; then
    printf '%s\n' "$PHASEPRED_DATA_ROOT"
  elif [[ -n "${XDG_DATA_HOME:-}" ]]; then
    printf '%s\n' "$XDG_DATA_HOME/phasepred"
  elif [[ -n "${HOME:-}" ]]; then
    printf '%s\n' "$HOME/.local/share/phasepred"
  else
    return 1
  fi
}

banner() {
  cat >&2 <<EOF
ESpritz is licensed under the Tosatto lab academic license v1.1, which forbids
redistribution. PhaSePred ships no ESpritz entity; you must obtain the archive
yourself.

  1. Request access from the Tosatto lab (contact details in the LICENSE file
     inside the archive). Their website is https://protein.bio.unipd.it/espritz/.
  2. Stage the downloaded archive as:
       \$PHASEPRED_VENDOR_ARCHIVE_DIR/$archive_name
     then rerun this installer. (A legacy copy at
     $pkg_dir/$archive_name is also accepted read-only.)
  3. Make sure perl is installed and on PATH.
EOF
}

check_available() {
  if [[ -n "${PHASEPRED_ESPRITZ_DIR:-}" ]] \
     && [[ -f "${PHASEPRED_ESPRITZ_DIR}/$entry_rel" ]]; then
    echo "ESpritz: ok (${PHASEPRED_ESPRITZ_DIR}/$entry_rel)"
    return 0
  fi
  if ! root="$(user_root)"; then
    echo "ESpritz: no data root (set HOME or PHASEPRED_DATA_ROOT)" >&2
    return 1
  fi
  if [[ -f "$root/espritz/$entry_rel" ]]; then
    echo "ESpritz: ok ($root/espritz/$entry_rel)"
    return 0
  fi
  if [[ -f "$pkg_dir/$entry_rel" ]]; then
    echo "ESpritz: ok ($pkg_dir/$entry_rel, vendored checkout)"
    return 0
  fi
  echo "ESpritz: MISSING ($entry_rel)" >&2
  echo "  Stage $archive_name in \$PHASEPRED_VENDOR_ARCHIVE_DIR and run" >&2
  echo "  bash $pkg_dir/install.sh" >&2
  return 1
}

if [[ "$mode" == "check" ]]; then
  check_available
  exit $?
fi

if ! root="$(user_root)"; then
  echo "ESpritz: no writable data root: set HOME or PHASEPRED_DATA_ROOT." >&2
  exit 1
fi
target="$root/espritz"

# Idempotency: an extracted user-data tree is never re-extracted.
if [[ -f "$target/$entry_rel" ]]; then
  echo "ESpritz: already installed ($target/$entry_rel) — skipping extraction."
  exec "$target/run" --check
fi

# docs/INSTALL.md §7: --offline is verify-only and must agree with --check,
# which is satisfied by the PHASEPRED_ESPRITZ_DIR override or the vendored
# checkout copy. Either exits 0 and creates nothing. Reverse proof: restoring
# the unconditional offline failure makes
# test_offline_agrees_with_check_for_every_installer fail.
if [[ "$offline" -eq 1 ]]; then
  if [[ -n "${PHASEPRED_ESPRITZ_DIR:-}" ]] \
     && [[ -f "${PHASEPRED_ESPRITZ_DIR}/$entry_rel" ]]; then
    echo "ESpritz: ok (${PHASEPRED_ESPRITZ_DIR}/$entry_rel) — nothing to download."
    exit 0
  fi
  if [[ -f "$pkg_dir/$entry_rel" ]]; then
    echo "ESpritz: ok ($pkg_dir/$entry_rel, vendored checkout) — nothing to download."
    exit 0
  fi
  if [[ ! -f "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}/$archive_name" ]]; then
    echo "ESpritz: --offline and no archive staged; cannot satisfy." >&2
    echo "  Set PHASEPRED_VENDOR_ARCHIVE_DIR to a directory holding $archive_name." >&2
    exit 1
  fi
fi

banner

# The only non-interactive source is the staged vendor archive. A checkout-local
# espritz.zip is read-only manual staging, accepted only after the interactive
# prompt: this is what makes the vendor branch load-bearing for `--offline` and
# `< /dev/null` (reverse proof T4-G4).
archive=""
if [[ -n "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}" ]] \
   && [[ -f "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name" ]]; then
  archive="${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name"
fi

if [[ -z "$archive" ]]; then
  echo "Press Enter once the archive is staged, or Ctrl+C to abort." >&2
  read -r _
  if [[ -n "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}" ]] \
     && [[ -f "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name" ]]; then
    archive="${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name"
  elif [[ -f "$pkg_dir/$archive_name" ]]; then
    archive="$pkg_dir/$archive_name"
  elif [[ -f "$target/$archive_name" ]]; then
    archive="$target/$archive_name"
  else
    echo "ESpritz: archive not found after prompt. Aborting." >&2
    exit 1
  fi
fi

if ! command -v unzip >/dev/null 2>&1; then
  echo "ESpritz: unzip not found on PATH." >&2
  exit 1
fi

if ! mkdir -p "$target" 2>/dev/null; then
  echo "ESpritz: cannot create $target (read-only home?)." >&2
  echo "  Set PHASEPRED_DATA_ROOT to a writable directory and retry." >&2
  exit 1
fi

echo "Extracting $archive into $target ..."
unzip -q -o "$archive" -d "$target"

if [[ ! -f "$target/$entry_rel" ]]; then
  echo "Unexpected archive layout — $entry_rel not found." >&2
  exit 1
fi

chmod +x "$target/espritz/espritz.pl"
chmod +x "$target/espritz/bin/"disbin* 2>/dev/null || true
chmod +x "$target/espritz/bin/blastredo" 2>/dev/null || true

cp "$pkg_dir/run" "$target/run"
chmod +x "$target/run"

echo "Done. Verify with:"
echo "  bash $pkg_dir/install.sh --check"
