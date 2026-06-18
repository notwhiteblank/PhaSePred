from __future__ import annotations

from pathlib import Path

from phasepred.legacy_features import (
    compute_lcr_features,
    compute_plaac_features,
    parse_plaac_output,
    parse_pscore_output,
)


def test_parse_pscore_output(tmp_path: Path) -> None:
    path = tmp_path / "pscore.tsv"
    path.write_text("PScore:                7.00                  >P1\n", encoding="utf-8")

    features = parse_pscore_output(path)

    assert features[0].accession == "P1"
    assert features[0].value == 7.0


def test_parse_plaac_output_uses_nllr() -> None:
    text = "SEQid\tMW\tNLLR\tPROTlen\nP1\t10\t0.217\t100\n"

    features = parse_plaac_output(text)

    assert features[0].accession == "P1"
    assert features[0].value == 0.217


def test_compute_lcr_features_from_seg() -> None:
    rows = [{"UniprotEntry": "P1", "sequence": "M" + "S" * 40}]

    features = compute_lcr_features(rows)

    assert features[0].accession == "P1"
    assert features[0].value > 0.5


def test_compute_plaac_features_smoke() -> None:
    rows = [{"UniprotEntry": "P1", "sequence": "Q" * 90}]

    features = compute_plaac_features(rows)

    assert features[0].accession == "P1"
    assert isinstance(features[0].value, float)
