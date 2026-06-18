from __future__ import annotations

import math

import pytest

from phasepred.catgranule_v1 import (
    TABLE_S4_SCALE_COLUMNS,
    TABLE_S4_SCALES,
    CatGranuleInputError,
    CatGranuleNormalization,
    CatGranuleScale,
    centered_window,
    granule_strength,
    normalize_score,
    profile,
    raw_score,
    residue_score,
    residue_scores,
)


def test_table_s4_scales_cover_standard_amino_acids() -> None:
    assert set(TABLE_S4_SCALES) == set("ACDEFGHIKLMNPQRSTVWY")
    assert TABLE_S4_SCALE_COLUMNS == ("rc", "rn", "dc", "dn", "pfg", "prg")


def test_table_s4_scales_match_supplement_examples() -> None:
    assert TABLE_S4_SCALES["A"] == CatGranuleScale(
        rc=0.45, rn=0.46, dc=0.08, dn=0.03, pfg=0.00, prg=0.00
    )
    assert TABLE_S4_SCALES["G"] == CatGranuleScale(
        rc=0.92, rn=0.52, dc=0.63, dn=1.00, pfg=1.00, prg=1.00
    )
    assert TABLE_S4_SCALES["R"] == CatGranuleScale(
        rc=0.89, rn=0.72, dc=0.17, dn=0.05, pfg=0.00, prg=1.00
    )
    assert TABLE_S4_SCALES["Y"] == CatGranuleScale(
        rc=1.00, rn=0.54, dc=0.45, dn=0.17, pfg=0.00, prg=0.00
    )


def test_centered_window_uses_truncated_heptapeptide_at_boundaries() -> None:
    sequence = "ACDEFGHIK"

    assert centered_window(sequence, 0) == "ACDE"
    assert centered_window(sequence, 4) == "CDEFGHI"
    assert centered_window(sequence, 8) == "GHIK"


def test_residue_score_averages_table_s4_scales_over_centered_window() -> None:
    # For ACG at index 1, the truncated window is ACG. The expected value is
    # the paper weights dotted against the mean Table S4 propensities.
    score = residue_score("ACG", 1)

    mean_rc = (0.45 + 0.00 + 0.92) / 3
    mean_rn = (0.46 + 0.00 + 0.52) / 3
    mean_dc = (0.08 + 0.24 + 0.63) / 3
    mean_dn = (0.03 + 0.01 + 1.00) / 3
    mean_prg = (0.00 + 0.00 + 1.00) / 3
    mean_pfg = (0.00 + 0.00 + 1.00) / 3
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
    sequence = "ACG"

    expected = sum(residue_scores(sequence)) / len(sequence) + 0.25 * math.log(len(sequence))

    assert math.isclose(raw_score(sequence), expected, rel_tol=1e-12)


def test_raw_score_can_omit_length_term_for_validation_experiments() -> None:
    sequence = "ACG"

    expected = sum(residue_scores(sequence)) / 3

    assert math.isclose(raw_score(sequence, include_length_term=False), expected)


def test_normalize_score_uses_external_mean_and_standard_deviation() -> None:
    stats = CatGranuleNormalization(mean=10.0, standard_deviation=2.0)

    assert normalize_score(13.0, stats) == 1.5


def test_normalization_rejects_non_positive_standard_deviation() -> None:
    with pytest.raises(CatGranuleInputError):
        CatGranuleNormalization(mean=0.0, standard_deviation=0.0)


def test_profile_averages_residue_scores_in_centered_window() -> None:
    sequence = "ACDEFG"
    scores = residue_scores(sequence)

    observed = profile(sequence, window_size=3)

    expected = [
        sum(scores[0:2]) / 2,
        sum(scores[0:3]) / 3,
        sum(scores[1:4]) / 3,
        sum(scores[2:5]) / 3,
        sum(scores[3:6]) / 3,
        sum(scores[4:6]) / 2,
    ]
    assert observed == pytest.approx(expected)


def test_profile_can_return_normalized_values() -> None:
    sequence = "ACDEFG"
    raw_profile = profile(sequence, window_size=3)
    stats = CatGranuleNormalization(mean=raw_profile[0], standard_deviation=2.0)

    observed = profile(sequence, window_size=3, normalization=stats)

    assert observed[0] == 0.0
    assert observed[1] == pytest.approx((raw_profile[1] - raw_profile[0]) / 2.0)


def test_granule_strength_averages_positive_normalized_profile_values() -> None:
    values = [-1.0, 0.0, 0.5, 1.5]

    assert granule_strength(values) == 1.0


def test_granule_strength_is_zero_when_profile_has_no_positive_values() -> None:
    assert granule_strength([-1.0, 0.0]) == 0.0


def test_scoring_rejects_empty_or_unsupported_sequences() -> None:
    with pytest.raises(CatGranuleInputError):
        raw_score("")

    with pytest.raises(CatGranuleInputError):
        raw_score("ACX")
