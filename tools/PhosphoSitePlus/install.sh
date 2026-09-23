#!/usr/bin/env bash
# PhosphoSitePlus installer — fetches the registration-gated dataset into the
# user data directory (idempotent, non-interactive).
#
#   bash install.sh             download (or reuse a vendor archive) and verify
#   bash install.sh --check     verify only; exit 0 present / 1 missing
#   bash install.sh --offline   never touch the network; exit 1 if unsatisfied
#
# Integrity policy (Ruling 3): PhosphoSitePlus is a rolling database, so its
# sha256 is a *known-version marker*. A mismatch prints a warning and records
# the actual digest, but installation still succeeds. This is the deliberate
# opposite of tools/PScore/install.sh, which hard-fails on a mismatch.
#
# Nothing is ever written into the checkout: the dataset and the install record
# land under the user data directory ($PHASEPRED_DATA_ROOT, else
# $XDG_DATA_HOME/phasepred, else ~/.local/share/phasepred).
set -euo pipefail

pkg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$pkg_dir/manifest.toml"
data_name="Phosphorylation_site_dataset.gz"

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

if [[ "$mode" == "check" ]]; then
  exec "$pkg_dir/run" --check
fi

if ! root="$(user_root)"; then
  echo "PhosphoSitePlus: no writable data root: set HOME or PHASEPRED_DATA_ROOT." >&2
  exit 1
fi
target_dir="$root/phosphositeplus"
target="$target_dir/$data_name"
url="$(manifest_get download_url)"
expected="$(manifest_get sha256)"
policy="$(manifest_get sha256_policy)"
declared_marker="$(manifest_get version_marker)"

if [[ -f "$target" ]]; then
  echo "PhosphoSitePlus: already installed ($target) — skipping download."
  exec "$pkg_dir/run" --check
fi

# docs/INSTALL.md §7: --offline is verify-only and must agree with --check.
# `run --check` is the exact predicate `install.sh --check` uses (it is the
# same exec), so a satisfied component — including the checkout's own
# data/raw/external copy — exits 0 before any mkdir. Reverse proof: restoring
# the unconditional offline failure makes
# test_offline_agrees_with_check_for_every_installer fail.
if [[ "$offline" -eq 1 ]]; then
  if "$pkg_dir/run" --check >/dev/null 2>&1; then
    echo "PhosphoSitePlus: already available — nothing to download."
    exec "$pkg_dir/run" --check
  fi
  if [[ ! -f "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}/$data_name" ]]; then
    echo "PhosphoSitePlus: --offline and no vendor archive; cannot satisfy." >&2
    echo "  Place $data_name in $target_dir or set PHASEPRED_VENDOR_ARCHIVE_DIR." >&2
    exit 1
  fi
fi

if ! mkdir -p "$target_dir" 2>/dev/null; then
  echo "PhosphoSitePlus: cannot create $target_dir (read-only home?)." >&2
  echo "  Set PHASEPRED_DATA_ROOT to a writable directory and retry." >&2
  exit 1
fi

staging="$target_dir/.$data_name.part"
rm -f "$staging"
source_desc=""
vendor_hit=0
if [[ -n "${PHASEPRED_VENDOR_ARCHIVE_DIR:-}" ]] \
   && [[ -f "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$data_name" ]]; then
  echo "PhosphoSitePlus: using vendor archive ${PHASEPRED_VENDOR_ARCHIVE_DIR}/$data_name"
  cp "${PHASEPRED_VENDOR_ARCHIVE_DIR}/$data_name" "$staging"
  source_desc="vendor:${PHASEPRED_VENDOR_ARCHIVE_DIR}/$data_name"
  vendor_hit=1
else
  echo "PhosphoSitePlus: downloading $url"
  if ! curl -fL --connect-timeout 15 --max-time 600 -o "$staging" "$url"; then
    rm -f "$staging"
    echo "PhosphoSitePlus: download failed (network or upstream error)." >&2
    echo "  Re-run, use --offline, or set PHASEPRED_VENDOR_ARCHIVE_DIR." >&2
    exit 1
  fi
  source_desc="$url"
fi

actual="$(sha256_of "$staging")"
mismatch=0
if [[ -n "$expected" && "$actual" != "$expected" ]]; then
  mismatch=1
  if [[ "$policy" == "hard" ]]; then
    rm -f "$staging"
    echo "PhosphoSitePlus: sha256 mismatch under hard policy — refusing to install." >&2
    echo "  expected $expected" >&2
    echo "  actual   $actual" >&2
    echo "  Nothing was written; fix the pin or the source." >&2
    exit 1
  fi
  echo "PhosphoSitePlus: WARNING sha256 mismatch (published policy: $policy)." >&2
  echo "  expected $expected" >&2
  echo "  actual   $actual" >&2
  echo "  Upstream has rolled; recording the actual digest and continuing." >&2
fi

mv -f "$staging" "$target"

marker=""
if command -v gzip >/dev/null 2>&1; then
  marker="$(gzip -dc "$target" 2>/dev/null | head -n1 || true)"
fi
if [[ -z "$marker" ]]; then
  marker="$declared_marker"
fi

{
  echo "component: PhosphoSitePlus"
  echo "installed_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "source: $source_desc"
  echo "path: $target"
  echo "sha256_expected: $expected"
  echo "sha256_actual: $actual"
  echo "sha256_policy: $policy"
  echo "sha256_match: $((1 - mismatch))"
  echo "version_marker: $marker"
} > "$target_dir/install-record.txt"

if [[ "$vendor_hit" -eq 1 ]]; then
  echo "PhosphoSitePlus: installed from vendor archive."
fi
echo "Done. Verify with: bash $pkg_dir/install.sh --check"
exec "$pkg_dir/run" --check
