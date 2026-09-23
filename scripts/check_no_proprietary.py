#!/usr/bin/env python3
"""Pre-commit guard: block non-redistributable and oversized files from git.

Scans ``git ls-files`` in the repository root and fails (exit 1) when a
tracked path matches a non-redistributable pattern (per
``docs/TOOL_LICENSES.md``) or exceeds a size threshold.

The patterns mirror ``.gitignore`` entries but also catch accidental
``git add -f`` and directory-rename mistakes (e.g. a re-introduced
``tools/per-tool/`` legacy layout).

Deliberately large, redistributable tracked data may be exempted from the size
cap via ``SIZE_ALLOWLIST`` (currently the ``data/processed/`` training tables).
The blocklist is checked first and always wins, so a non-redistributable path is
never exempted by size.

Usage (run from the repository root):

    python scripts/check_no_proprietary.py            # default 5 MB cap
    python scripts/check_no_proprietary.py --large-mb 10
    python scripts/check_no_proprietary.py --quiet    # only affected paths
"""

from __future__ import annotations

import argparse
import fnmatch
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# Non-redistributable / not-to-be-shipped tracked paths (POSIX-relative to
# the repo root, fnmatch style). Sources: docs/TOOL_LICENSES.md §1.
BLOCKLIST = [
    # Legacy three-way layout must stay gone (S3 normalized into contracts).
    "tools/per-tool/*",
    "tools/wrappers/*",
    "tools/install/*",
    # ESpritz (Tosatto academic license — no redistribution).
    "tools/ESpritz/espritz/*",
    "tools/ESpritz/*.zip",
    # PScore data (eLife CC-BY code ships; DBS data stays out of git).
    "tools/PScore/*.tgz",
    "tools/PScore/SourceCodeS2/DBS/*",
    # PLAAC upstream archive (MIT code ships; big zip stays out).
    "tools/PLAAC/*.zip",
    # IUPred3 (ELTE academic license — archive/entity not redistributable).
    "tools/_archive/IUPred3/iupred3/*",
    "tools/_archive/IUPred3/*.tar.gz",
    "tools/_archive/*/iupred3.tar.gz",
    # catGRANULE v2 (experiment archive).
    "tools/_archive/catGRANULE_v2/*.tar.gz",
    "tools/_archive/catGRANULE_v2/catGRANULE2.0-1.0.0/*",
    # Isolated conda environments (DeepCoil / catGRANULE2).
    ".external_envs/*",
    # PhosphoSitePlus registration-only dataset.
    "data/raw/external/phosphositeplus/*",
    # E4 data contracts + user-data roots. If PHASEPRED_DATA_ROOT points inside
    # the checkout these hold non-redistributable downloads (Ruling 6). The
    # bundled DeepPhase TSV/LICENSE/NOTICE (MIT) live under
    # src/phasepred/data/deepphase/ and are deliberately NOT matched.
    "tools/PhosphoSitePlus/*.gz",
    "phosphositeplus/*",
    "pscore/*",
    "espritz/*",
    "envs/*",
]

# Tracked paths exempt from the size cap. These are deliberate, documented
# artifacts derived from the paper supplementary datasets (E1): the four
# `chen2022_s2s3_*.tsv` training tables plus `MANIFEST.tsv` and `README.md`. The
# largest is 12.8 MiB. The cap exists to catch accidental blobs, not these. The
# blocklist above is checked first and still wins, so a non-redistributable path
# is never exempted by size.
SIZE_ALLOWLIST = [
    "data/processed/chen2022_s2s3_*.tsv",
    "data/processed/MANIFEST.tsv",
    "data/processed/README.md",
]


def tracked_paths() -> list[str]:
    proc = subprocess.run(
        ["git", "ls-files", "-z"], cwd=str(REPO_ROOT), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"git ls-files failed: {proc.stderr}")
    return [p for p in proc.stdout.split("\0") if p]


def _matches_at_any_depth(path: str, pattern: str) -> bool:
    """Match ``pattern`` against the path or any of its suffix segment runs.

    A root-anchored pattern like ``pscore/*`` must also catch a user data root
    nested deeper inside the checkout (e.g. ``vendor/pscore/...``), so every
    descending tail of the path is tried. Suffix tails are whole path segments,
    so ``pscore/*`` matches ``a/b/pscore/x`` but not ``a/b/myscore/x``.
    """
    parts = path.split("/")
    for i in range(len(parts)):
        tail = "/".join(parts[i:])
        if fnmatch.fnmatch(tail, pattern) or fnmatch.fnmatch(
            tail, pattern.rstrip("/") + "*"
        ):
            return True
    return False


def matches_blocklist(path: str) -> bool:
    for pattern in BLOCKLIST:
        if _matches_at_any_depth(path, pattern):
            return True
    return False


def matches_size_allowlist(path: str) -> bool:
    for pattern in SIZE_ALLOWLIST:
        if fnmatch.fnmatch(path, pattern) or fnmatch.fnmatch(path, pattern.rstrip("/") + "*"):
            return True
    return False


def size_mb(path: str) -> float:
    try:
        return (REPO_ROOT / path).stat().st_size / (1024 * 1024)
    except OSError:
        return 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large-mb", type=float, default=5.0, help="size cap in MiB")
    parser.add_argument("--quiet", action="store_true", help="print only offending paths")
    args = parser.parse_args()

    violations: list[str] = []
    oversized: list[str] = []
    paths = tracked_paths()
    for path in paths:
        if matches_blocklist(path):
            violations.append(path)
        elif not matches_size_allowlist(path) and size_mb(path) > args.large_mb:
            oversized.append(f"{path}  ({size_mb(path):.1f} MiB)")

    clean = not violations and not oversized
    if args.quiet:
        for path in violations + oversized:
            print(path)
    else:
        if violations:
            print("FATAL: tracked paths match non-redistributable patterns:")
            for path in violations:
                print(f"  - {path}")
        if oversized:
            print(f"FATAL: tracked files exceed --large-mb={args.large_mb:g}:")
            for path in oversized:
                print(f"  - {path}")
        print(f"\nscanned {len(paths)} tracked files: "
              f"{len(violations)} blocked pattern(s), {len(oversized)} oversized")
        if clean:
            print("OK — no non-redistributable or oversized files tracked.")

    return 0 if clean else 1


if __name__ == "__main__":
    sys.exit(main())
