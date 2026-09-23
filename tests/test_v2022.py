"""v2022 feature-definition tests (docs/FEATURES.md).

Each deviation that S2 fixes gets a unit test here with expectations
anchored to the S1 web fixtures (Q08211 / P35637) where applicable.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phasepred.features import (
    CATGRANULE_RESIDUE_OFFSET,
    CATGRANULE_RESIDUE_WINDOW,
    DEEPCOIL_THRESHOLD,
    binary_deepcoil,
    catgranule_profile_alignment,
    pscore_alignment,
)
from phasepred.tools import _binarize_deepcoil_frame, _parse_espritz_states


def _q08211_fasta() -> str:
    return Path(
        __file__, "..", "..", "tests", "fixtures", "sequences", "Q08211.fasta"
    ).resolve().read_text(encoding="utf-8")


def _q08211_sequence() -> str:
    return "".join(
        line
        for line in _q08211_fasta().splitlines()
        if not line.startswith(">")
    )


# ---------------------------------------------------------------------------
# DeepCoil v2022 definition: binary threshold 0.82 (FEATURES.md §8)
# ---------------------------------------------------------------------------


def test_deepcoil_binary_below_threshold() -> None:
    raw_cc = [0.0, 0.0107, 0.8, 0.819]
    assert binary_deepcoil(raw_cc) == 0.0


def test_deepcoil_binary_at_threshold() -> None:
    assert binary_deepcoil([0.82]) == 1.0
    assert binary_deepcoil([0.819, 0.8200001]) == 1.0
    # Ruling 9(c): pin the >= operator on the exact literals the archive and the
    # local O75146 run straddle (archive 0.8199 -> 0, local 0.82 -> 1,
    # 0.8201 -> 1). Synthetic: no protein, environment or archive needed.
    assert binary_deepcoil([0.8199]) == 0.0
    assert binary_deepcoil([0.82]) == 1.0
    assert binary_deepcoil([0.8201]) == 1.0
    assert DEEPCOIL_THRESHOLD == 0.82


def test_deepcoil_binary_frame_uses_max_raw_cc() -> None:
    frame = pd.DataFrame({"raw_cc": [0.1, 0.9, 0.2]})
    assert _binarize_deepcoil_frame(frame) == 1.0

    frame = pd.DataFrame({"raw_cc": [0.1, 0.81, 0.2]})
    assert _binarize_deepcoil_frame(frame) == 0.0


def test_deepcoil_q08211_anchor() -> None:
    """Q08211 web deepcoil.single==0.0; local raw_cc max==0.0107 (<0.82)."""
    assert binary_deepcoil([0.0107]) == 0.0
    assert DEEPCOIL_THRESHOLD == 0.82


# ---------------------------------------------------------------------------
# ESpritz IDR parsing fix (FEATURES.md §3): banner lines must not count
# ---------------------------------------------------------------------------


def test_parse_espritz_states_ignores_banner_and_license_lines() -> None:
    text = """\
################################################################################################

Licensed to: PhD Kaiqiang You (Academic) — this license is for non-commercial use only.
 Please read our LICENSE file for details.

model : D
thres : 0.5072
/tmp/esp_out/Q08211.fasta
Finished executing DISPROT disorder no psi-blast and threshold=0.5072

Time for predictions to finish = 0 hrs : 0 mins : 0 secs
******************************************************************************************************
O\t0.27802
O\t0.280725
D\t0.81
O\t0.21
"""
    path = Path("/tmp") / "q08211_states.espritz"
    path.write_text(text, encoding="utf-8")

    try:
        states = _parse_espritz_states(path)
    finally:
        path.unlink(missing_ok=True)

    assert states == ["O", "O", "D", "O"]
    idr = sum(s == "D" for s in states) / len(states)
    assert idr == 0.25


def test_espritz_idr_q08211_anchor_is_zero() -> None:
    """Q08211 web espritz.single==0.0 (all-O label); the old banner-counting
    parse produced 7.2e-05, the fixed two-column parse must give 0.0."""
    web_path = (
        Path(__file__).resolve().parents[1]
        / "tests" / "fixtures" / "sequences" / "Q08211.fasta"
    )
    assert web_path.exists()
    label_path = Path("/tmp") / "q08211_anchor.espritz"
    label_path.write_text("model : D\nthres : 0.5072\n" + "O\t0.1\n" * 1270, encoding="utf-8")
    try:
        states = _parse_espritz_states(label_path)
    finally:
        label_path.unlink(missing_ok=True)
    assert len(states) == 1270
    assert sum(s == "D" for s in states) == 0


# ---------------------------------------------------------------------------
# Residue-level alignment utilities (FEATURES.md §5/§7)
# ---------------------------------------------------------------------------


def test_catgranule_profile_alignment_q08211() -> None:
    positions = catgranule_profile_alignment(1270)
    assert len(positions) == 1220  # L - 50
    assert positions.start == CATGRANULE_RESIDUE_OFFSET == 24
    assert positions[0] == 24
    assert positions[-1] == 1270 - 27  # last interior center (L-1(0-idx) - 25 - 1)
    assert CATGRANULE_RESIDUE_WINDOW == 51


def test_catgranule_profile_alignment_short_sequence() -> None:
    assert list(catgranule_profile_alignment(10)) == []


def test_pscore_alignment_q08211_suffix() -> None:
    # web pscore residue length 1267 = L - 3 (suffix aligned)
    positions = pscore_alignment(1270, 1267)
    assert len(positions) == 1267
    assert positions.start == 3
    assert positions[-1] == 1269


def test_pscore_alignment_oversized_array() -> None:
    assert list(pscore_alignment(10, 11)) == []


# ---------------------------------------------------------------------------
# Native v2022 anchors (QD08211 web values vs local computation)
# ---------------------------------------------------------------------------


def test_native_feature_anchors_q08211() -> None:
    from phasepred.features import compute_native_features

    features = compute_native_features(_q08211_sequence())

    assert features["length"] == 1270
    assert features["Hydropathy"] == pytest.approx(0.462353455818, abs=1e-9)
    assert features["FCR"] == pytest.approx(0.225196850394, abs=1e-9)