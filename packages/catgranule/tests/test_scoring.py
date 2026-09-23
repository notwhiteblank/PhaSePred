from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest
from catgranule import score_sequence
from catgranule.scoring import (
    CatGranuleInputError,
    centered_window,
    granule_strength,
    normalize,
    raw_score,
    residue_profile,
    residue_score,
    residue_scores,
)
from catgranule.weights import (
    SCALE_COLUMNS,
    CatGranuleWeightsError,
    load_weights,
)

FIXTURE_FASTA = Path(__file__).parent / "fixtures" / "Q08211.fasta"
FIXTURE_P35637 = Path(__file__).parent / "fixtures" / "P35637.fasta"

# Web-archive anchors (tests/fixtures/golden/web_*.json) and the holdout
# |delta| distribution they pass through (distilled: med 0.0177, q90 0.057).
WEB_Q08211_CG = 1.43649
WEB_P35637_CG = 5.750
DISTILLED_TOL = 0.05

ROUTES = ("distilled", "paper", "legacy-unaudited")


def _q08211() -> str:
    return "".join(
        line
        for line in FIXTURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )


def _p35637() -> str:
    return "".join(
        line
        for line in FIXTURE_P35637.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )


@pytest.mark.parametrize("route", ROUTES)
def test_scales_cover_standard_amino_acids_all_routes(route: str) -> None:
    weights = load_weights(route=route)
    assert set(weights.scales) == set("ACDEFGHIKLMNPQRSTVWY")
    assert SCALE_COLUMNS == ("rc", "rn", "dc", "dn", "pfg", "prg")


def test_default_weights_are_distilled() -> None:
    weights = load_weights()
    assert weights.provenance.route == "distilled"
    assert weights.mode == "distilled"
    assert weights.model_path is not None
    assert weights.model_path.exists()


@pytest.mark.parametrize("route", ROUTES)
def test_scales_match_mmc1_table_s4_examples(route: str) -> None:
    weights = load_weights(route=route)
    assert weights.scales["A"] == {
        "rc": 0.45, "rn": 0.46, "dc": 0.08, "dn": 0.03, "pfg": 0.0, "prg": 0.0,
    }
    assert weights.scales["G"] == {
        "rc": 0.92, "rn": 0.52, "dc": 0.63, "dn": 1.00, "pfg": 1.0, "prg": 1.0,
    }
    assert weights.scales["R"] == {
        "rc": 0.89, "rn": 0.72, "dc": 0.17, "dn": 0.05, "pfg": 0.0, "prg": 1.0,
    }


def test_centered_window_uses_truncated_heptapeptide_at_boundaries() -> None:
    weights = load_weights()
    sequence = "ACDEFGHIK"

    assert centered_window(sequence, 0, radius=3, scales=weights.scales) == "ACDE"
    assert centered_window(sequence, 4, radius=3, scales=weights.scales) == "CDEFGHI"
    assert centered_window(sequence, 8, radius=3, scales=weights.scales) == "GHIK"


def test_residue_score_averages_scales_over_centered_window() -> None:
    # For ACG at index 1, the truncated window is ACG. The expected value is
    # the paper weights dotted against the mean residue propensities.
    weights = load_weights()
    score = residue_score("ACG", 1, weights)

    mean_rc = (0.45 + 0.0 + 0.92) / 3
    mean_rn = (0.46 + 0.0 + 0.52) / 3
    mean_dc = (0.08 + 0.24 + 0.63) / 3
    mean_dn = (0.03 + 0.01 + 1.00) / 3
    mean_prg = (0.0 + 0.0 + 1.00) / 3
    mean_pfg = (0.0 + 0.0 + 1.00) / 3
    expected = (
        0.48 * mean_rc
        + 7.24 * mean_rn
        + 0.26 * mean_dc
        + 11.54 * mean_dn
        + 1.98 * mean_prg
        + 1.42 * mean_pfg
    )

    assert math.isclose(score, expected, rel_tol=1e-12)


def test_raw_score_is_mean_residue_score_plus_log_length_term() -> None:
    weights = load_weights()
    sequence = "ACG"

    expected = sum(residue_scores(sequence, weights)) / len(sequence) + 0.25 * math.log(
        len(sequence)
    )

    assert math.isclose(raw_score(sequence, weights), expected, rel_tol=1e-12)


def test_raw_score_can_omit_length_term() -> None:
    weights = load_weights()
    sequence = "ACG"

    expected = sum(residue_scores(sequence, weights)) / 3

    assert math.isclose(raw_score(sequence, weights, include_length_term=False), expected)


def test_normalize_uses_artifact_constants() -> None:
    weights = load_weights()

    assert normalize(weights.normalization_mean, weights) == pytest.approx(0.0)
    assert normalize(
        weights.normalization_mean + weights.normalization_std, weights
    ) == pytest.approx(1.0)


def test_residue_profile_covers_interior_only() -> None:
    weights = load_weights()

    profile = residue_profile("ACDEFGHIK", weights)

    # 9 residues, 51-window interior impossible -> no covered positions
    assert profile == []


def test_residue_profile_length_is_l_minus_window_plus_one() -> None:
    weights = load_weights()
    sequence = "A" * 120

    profile = residue_profile(sequence, weights)

    assert len(profile) == 120 - (weights.residue_window_size - 1)


def test_granule_strength_averages_positive_normalized_profile_values() -> None:
    assert granule_strength([-1.0, 0.0, 0.5, 1.5]) == 1.0


def test_granule_strength_is_zero_when_profile_has_no_positive_values() -> None:
    assert granule_strength([-1.0, 0.0]) == 0.0


def test_score_sequence_rejects_empty_or_unsupported_sequences() -> None:
    with pytest.raises(CatGranuleInputError):
        score_sequence("")

    with pytest.raises(CatGranuleInputError):
        score_sequence("ACX")


def test_q08211_distilled_anchor() -> None:
    """Q08211 anchor — distilled (default) weights must approximate the web-archive
    single (1.43649) within the holdout |delta| distribution (med 0.0177, q90 0.057)."""
    scored = score_sequence(_q08211())
    assert abs(scored["single"] - WEB_Q08211_CG) <= DISTILLED_TOL
    assert len(scored["residue"]) == 1220  # L - 50


def test_p35637_distilled_anchor() -> None:
    """P35637 anchor — distilled (default) weights vs web-archive 5.750."""
    scored = score_sequence(_p35637())
    assert abs(scored["single"] - WEB_P35637_CG) <= DISTILLED_TOL
    assert len(scored["residue"]) == 476  # L - 50


def test_legacy_anchor_continuity() -> None:
    """legacy-unaudited route preserves the S0/S1 local reconstruction value."""
    scored = score_sequence(_q08211(), route="legacy-unaudited")
    assert math.isclose(scored["single"], 1.5206734862994171, rel_tol=1e-9)
    assert len(scored["residue"]) == 1220  # L - 50


def test_paper_anchor_q08211() -> None:
    """paper route reproduces the S4 main-route normalized single for Q08211."""
    scored = score_sequence(_q08211(), route="paper")
    assert math.isclose(scored["single"], 1.517035, rel_tol=1e-5)


def test_distilled_residue_profile_is_aligned_to_single() -> None:
    """For the distilled route the residue profile is translated so its mean
    equals the surrogate single (per-protein affine alignment)."""
    scored = score_sequence(_q08211())
    profile = scored["residue"]
    assert profile
    mean_profile = sum(profile) / len(profile)
    assert mean_profile == pytest.approx(scored["single"], rel=1e-9)


def test_q08211_web_profile_correlation() -> None:
    """The 51-window interior profile must track the web archive residue
    array (Pearson ~0.999 per S1 evidence)."""
    import numpy as np

    repo_root = Path(__file__).resolve().parents[3]
    web_path = repo_root / "tests" / "fixtures" / "golden" / "web_Q08211.json"
    if not web_path.exists():
        override = Path(os.environ.get("PHASEPRED_DATA_ROOT", str(repo_root)))
        web_path = override / "tests" / "fixtures" / "golden" / "web_Q08211.json"
        if not web_path.exists():
            pytest.skip("web_Q08211.json fixture not available outside the PhaSePred repo")
    web = json.loads(web_path.read_text(encoding="utf-8"))
    local = score_sequence(_q08211())["residue"]
    a = np.array(local)
    b = np.array(web["data"]["tools"]["catgranule"]["residue"])
    assert len(a) == len(b) == 1220
    assert np.corrcoef(a, b)[0, 1] > 0.99


def test_weights_provenance_is_distilled_default() -> None:
    weights = load_weights()

    assert not weights.is_legacy_unaudited
    assert weights.is_distilled
    assert weights.provenance.route == "distilled"
    assert weights.provenance.gate is not None
    assert weights.provenance.source


def test_weights_paper_provenance() -> None:
    weights = load_weights(route="paper")
    assert weights.provenance.route == "paper"
    assert weights.provenance.gate is not None
    pos = (
        "alternate artifact; original S4 main route "
        "(below gate, see CATGRANULE_VALIDATION.md)"
    )
    assert weights.provenance.gate["position"] == pos


def test_weights_legacy_provenance() -> None:
    weights = load_weights(route="legacy-unaudited")
    assert weights.is_legacy_unaudited
    assert weights.provenance.gate is None


def test_weights_artifact_schema_validation() -> None:
    weights = load_weights()

    raw = weights.to_json()
    assert raw["provenance"]["route"] == "distilled"
    assert raw["mode"] == "distilled"
    assert raw["model_file"] == "distilled_xgb_model.ubj"
    assert set(raw["weights"]) == set(SCALE_COLUMNS) | {"length"}
    assert len(raw["scales"]) == 20


def test_weights_override_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from catgranule.weights import _read_artifact, resolve_artifact_path

    artifact = _read_artifact(resolve_artifact_path(route="paper"))
    artifact["provenance"]["route"] = "distilled"
    artifact["provenance"]["gate"] = {"rho": 0.99, "delta_median": 0.004, "holdout": "web-human"}
    overridden = tmp_path / "weights.json"
    overridden.write_text(json.dumps(artifact), encoding="utf-8")
    monkeypatch.setenv("CATGRANULE_WEIGHTS", str(overridden))

    weights = load_weights()

    assert weights.provenance.route == "distilled"
    assert weights.provenance.gate == {"rho": 0.99, "delta_median": 0.004, "holdout": "web-human"}


def test_weights_rejects_bad_provenance(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from catgranule.weights import _read_artifact, resolve_artifact_path

    artifact = _read_artifact(resolve_artifact_path(route="paper"))
    artifact["provenance"]["route"] = "mystery"
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(artifact), encoding="utf-8")
    monkeypatch.setenv("CATGRANULE_WEIGHTS", str(bad))

    with pytest.raises(CatGranuleWeightsError):
        load_weights()


def test_score_sequence_route_param_switches_artifacts() -> None:
    q = _q08211()
    distilled = score_sequence(q, route="distilled")
    paper = score_sequence(q, route="paper")
    legacy = score_sequence(q, route="legacy-unaudited")
    # surrogate (web-archive caliber), paper and legacy are all distinct values
    assert distilled["single"] != pytest.approx(paper["single"])
    assert paper["single"] != pytest.approx(legacy["single"])
    # formula routes agree on the raw (un-normalized) basis
    assert raw_score(q, load_weights(route="paper")) == pytest.approx(
        raw_score(q, load_weights(route="legacy-unaudited"))
    )


def test_score_batch_iterates_in_order() -> None:
    from catgranule import score_batch

    records = [
        {"accession": "P1", "sequence": "ACD"},
        {"accession": "P2", "sequence": "ACDEFGHIK"},
    ]
    scored = list(score_batch(records))
    assert [item["accession"] for item in scored] == ["P1", "P2"]
    assert set(scored[0]) == {"accession", "single", "residue", "granule_strength"}
    assert scored[0]["single"] == score_sequence("ACD")["single"]


def test_score_batch_accepts_a_generator() -> None:
    from catgranule import score_batch

    def gen():
        yield {"accession": "P1", "sequence": "ACD"}
        yield {"accession": "P2", "sequence": "ACDEFGHIK"}

    scored = list(score_batch(gen()))
    assert [item["accession"] for item in scored] == ["P1", "P2"]


def test_route_argument_matches_weights_overload() -> None:
    q = _q08211()
    via_route = score_sequence(q, route="paper")
    via_weights = score_sequence(q, weights=load_weights(route="paper"))
    assert via_route["single"] == via_weights["single"]


def test_distilled_feature_map_matches_s4_numeric() -> None:
    """The shipped feature map must be bit-identical to the training-time map."""
    import sys

    import numpy as np

    scripts_dir = Path(__file__).resolve().parents[1] / "scripts"
    sys.path.insert(0, str(scripts_dir))
    try:
        import _s4_numeric
        from catgranule._distill import (
            coefficients as shipped_coeffs,
        )
        from catgranule._distill import (
            featurize as shipped,
        )
        from catgranule._distill import (
            scale_matrix as shipped_matrix,
        )
        from catgranule.weights import load_weights

        w = load_weights(route="legacy-unaudited")
        sm = shipped_matrix(w)
        co = shipped_coeffs(w)
        probe = "MANEPQNNDSGLYFRILLQVEADEPTINAIVNSGSDISDAVQGISSNQNIISGNESRGG"[:60]
        for seq in (_q08211(), _p35637(), probe):
            a = shipped(seq, w)
            b = _s4_numeric.featurize_distill(seq, sm, co)
            assert a.shape == b.shape == (181,)
            assert np.max(np.abs(a - b)) < 1e-12, np.max(np.abs(a - b))
    finally:
        sys.path.pop(0)


def test_score_batch_raises_on_malformed_record() -> None:
    from catgranule import score_batch

    with pytest.raises(CatGranuleInputError):
        list(score_batch([{"id": "P1", "seq": "ACD"}]))
    with pytest.raises(CatGranuleInputError):
        list(score_batch([{"accession": "P1"}]))


def test_q08211_anchor_raw_value_stable() -> None:
    weights = load_weights(route="legacy-unaudited")
    raw = raw_score(_q08211(), weights)

    # raw Equation-2 value captured from S1 reconstruction
    target = weights.normalization_mean + weights.normalization_std * 1.5206734862994171
    assert math.isclose(raw, target, rel_tol=1e-6)