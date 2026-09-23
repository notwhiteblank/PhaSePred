#!/usr/bin/env bash
# PScore installer — fetches the eLife Source-code-2 archive, verifies it, and
# unpacks SourceCodeS2/ + DBS/ into the user data directory (idempotent).
#
#   bash install.sh             download (or reuse a vendor archive) + verify
#   bash install.sh --check     verify only; exit 0 available / 1 missing
#   bash install.sh --offline   never touch the network; exit 1 if unsatisfied
#
# Integrity policy (Ruling 3): the eLife tgz is a static 2019 archive, so its
# sha256 is *hard* — a mismatch deletes the partial download, removes any
# half-extracted tree, prints a clear error, and exits non-zero. This is the
# deliberate opposite of tools/PhosphoSitePlus/install.sh (known-version marker).
#
# Nothing is written into the checkout: the archive, SourceCodeS2/, DBS/, and
# the copied `run` all land under the user data directory
# ($PHASEPRED_DATA_ROOT, else $XDG_DATA_HOME/phasepred, else
# ~/.local/share/phasepred).
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$pkg_dir/manifest.toml"
archive_name="elife-31486-code2-v2.tgz"
predictor_rel="SourceCodeS2/elife_phase_separation_predictor.py"
dbs_rel="SourceCodeS2/DBS"

mode="install"
offline=0
for arg in "$@"; do
  case "$arg" in
    --check) mode="check" ;;
    --offline) offline=1 ;;
    -h|--help)
      sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Usage: bash $0 [--check] [--offline]" >&2
      exit 2
      ;;
  esac
done

manifest_get() {
  sed -n "s/^$1 *= *\"\\(.*\\)\"$/\\1/p" "$manifest" | head -n1
}

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

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

available_in() {
  local dir="$1"
  [[ -f "$dir/$predictor_rel" && -d "$dir/$dbs_rel" ]]
}

if ! root="$(user_root)"; then
  if [[ "$mode" == "check" ]]; then
    echo "PScore: no data root (set HOME or PHASEPRED_DATA_ROOT)" >&2
    exit 1
  fi
  echo "PScore: no writable data root: set HOME or PHASEPRED_DATA_ROOT." >&2
  exit 1
fi
target="$root/pscore"
url="$(manifest_get download_url)"
expected="$(manifest_get sha256)"
policy="$(manifest_get sha256_policy)"

if [[ "$mode" == "check" ]]; then
  if available_in "$target"; then
    echo "PScore: ok ($target)"
    exec "$pkg_dir/run" --check
  fi
  if available_in "$pkg_dir"; then
    echo "PScore: ok ($pkg_dir, vendored checkout)"
    exec "$pkg_dir/run" --check
  fi
  echo "PScore: MISSING ($predictor_rel + $dbs_rel)" >&2
  echo "  Install it with: bash $pkg_dir/install.sh" >&2
  echo "  User data target: $target" >&2
  exit 1
fi

# Idempotency: a complete user-data copy is never re-downloaded.
if available_in "$target"; then
  if [[ ! -f "$target/run" ]]; then
    cp "$pkg_dir/run" "$target/run"
    chmod +x "$target/run"
  fi
  echo "PScore: already installed ($target/SourceCodeS2) — skipping download."
  exec "$target/run" --check
fi

# docs/INSTALL.md §7: --offline is a verify-only no-op when the component is
# already satisfied, and it must not create anything. The vendored checkout copy
# counts as satisfied (exactly as --check does), so an offline run in a checkout
# that holds it exits 0 before any mkdir. Reverse proof: restoring the old order
# (mkdir, then an unconditional offline failure) makes
# test_offline_agrees_with_check_for_every_installer fail.
if [[ "$offline" -eq 1 ]]; then
  if available_in "$pkg_dir"; then
    echo "PScore: ok ($pkg_dir, vendored checkout) — nothing to download."
    exec "$pkg_dir/run" --check
  fi
  if [[ ! -f "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}/$archive_name" ]]; then
    echo "PScore: --offline and no vendor archive; cannot satisfy." >&2
    echo "  Place $archive_name in $target or set PHASEPRED_VENDOR_ARCHIVE_DIR." >&2
    exit 1
  fi
fi

if ! mkdir -p "$target" 2>/dev/null; then
  echo "PScore: cannot create $target (read-only home?)." >&2
  echo "  Set PHASEPRED_DATA_ROOT to a writable directory and retry." >&2
  exit 1
fi

staging="$target/.$archive_name.part"
rm -f "$staging"
vendor_hit=0
if [[ -n "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}" ]] \
   && [[ -f "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name" ]]; then
  echo "PScore: using vendor archive ${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name"
  cp "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$archive_name" "$staging"
  vendor_hit=1
else
  echo "PScore: downloading $url"
  if ! curl -fL --connect-timeout 15 --max-time 600 -o "$staging" "$url"; then
    rm -f "$staging"
    echo "PScore: download failed (network or upstream error)." >&2
    echo "  Re-run, use --offline, or set PHASEPRED_VENDOR_ARCHIVE_DIR." >&2
    exit 1
  fi
fi

actual="$(sha256_of "$staging")"
if [[ -n "$expected" && "$actual" != "$expected" ]]; then
  if [[ "$policy" == "known-version" ]]; then
    echo "PScore: WARNING sha256 mismatch (published policy: $policy)." >&2
    echo "  expected $expected" >&2
    echo "  actual   $actual" >&2
  else
    rm -f "$staging"
    rm -rf "$target/SourceCodeS2"
    echo "PScore: sha256 mismatch — refusing to install." >&2
    echo "  expected $expected" >&2
    echo "  actual   $actual" >&2
    echo "  The partial download and any half-extracted tree were removed." >&2
    exit 1
  fi
fi

mv -f "$staging" "$target/$archive_name"

extract_dir="$target/.extract.$$"
rm -rf "$extract_dir"
mkdir -p "$extract_dir"
if ! tar -xzf "$target/$archive_name" -C "$extract_dir"; then
  rm -rf "$extract_dir"
  echo "PScore: archive extraction failed; no partial tree was kept." >&2
  exit 1
fi
if [[ ! -f "$extract_dir/$predictor_rel" || ! -d "$extract_dir/$dbs_rel" ]]; then
  rm -rf "$extract_dir"
  echo "PScore: unexpected archive layout — $predictor_rel + $dbs_rel not found." >&2
  exit 1
fi

rm -rf "$target/SourceCodeS2"
mv "$extract_dir/SourceCodeS2" "$target/SourceCodeS2"
rm -rf "$extract_dir"

cp "$pkg_dir/run" "$target/run"
chmod +x "$target/run"

if [[ "$vendor_hit" -eq 1 ]]; then
  echo "PScore: installed from vendor archive."
fi
echo "Done. Verify with: bash $pkg_dir/install.sh --check"
exec "$target/run" --check
