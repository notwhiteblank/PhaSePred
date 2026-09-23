#!/usr/bin/env python
"""Product A re-training (paper protocol, v2022).

This module is the thin, canonical location wrapper around the authoritative
one-command retrain script ``scripts/train_phasepred.py`` (see
RETRAIN_PROTOCOL.md §6: the train script location stays in
``products/A_paper_split_recomputed/``; the implementation reads the committed
``data/processed/`` tables and follows the paper's default-XGB 5-fold x
5-round x 10 negative-set protocol). A bare ``runpy`` of the target matches the
wrapper's original behaviour: the target's default output directories are
``src/phasepred/data/models`` and ``products/A_paper_split_recomputed``.
"""

from __future__ import annotations

import runpy
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

runpy.run_path(str(REPO_ROOT / "scripts" / "train_phasepred.py"), run_name="__main__")
