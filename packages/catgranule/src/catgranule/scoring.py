"""catGRANULE scoring: Bolognesi et al. 2016, Cell Reports 16, equations 1-2.

Features for every residue use a centered heptapeptide window over the six
normalized physico-chemical scales (RC/RN/DC/DN/PRG/PFG; mmc1 Table S4).
The protein-level ``single`` score is the mean residue score plus a length
term, then optionally Z-normalized with the constants shipped in the
weights artifact. The residue-level profile is a 51-residue sliding-window
mean (radius 25) restricted to the fully-covered interior positions, which
is the length-``L-50`` form used by the PhaSePred web archive
(often written as "offset 24, length L-50").
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Iterator, Mapping

from catgranule.weights import SCALE_COLUMNS, CatGranuleWeights, load_weights

_AA_IDENTIFIER_RE = re.compile(r"^[A-Za-z]+$")


class CatGranuleInputError(ValueError):
    """Raised when a sequence cannot be scored."""


def centered_window(
    sequence: str, index: int, *, radius: int, scales: Mapping[str, object]
) -> str:
    normalized = normalize_sequence(sequence, scales)
    if index < 0 or index >= len(normalized):
        raise CatGranuleInputError(f"Sequence index out of range: {index}")
    start = max(0, index - radius)
    stop = min(len(normalized), index + radius + 1)
    return normalized[start:stop]


def residue_score(sequence: str, index: int, weights: CatGranuleWeights) -> float:
    window = centered_window(sequence, index, radius=weights.window_radius, scales=weights.scales)
    width = len(window)
    means = {
        column: sum(weights.scales[residue][column] for residue in window) / width
        for column in SCALE_COLUMNS
    }
    return sum(float(weights.weights[column]) * means[column] for column in SCALE_COLUMNS)


def residue_scores(sequence: str, weights: CatGranuleWeights) -> list[float]:
    normalized = normalize_sequence(sequence, weights.scales)
    return [residue_score(normalized, index, weights) for index in range(len(normalized))]


def raw_score(
    sequence: str, weights: CatGranuleWeights, *, include_length_term: bool = True
) -> float:
    """Equation 2: mean residue score plus a log-length term (Bolognesi 2016)."""
    normalized = normalize_sequence(sequence, weights.scales)
    scores = residue_scores(normalized, weights)
    score = sum(scores) / len(scores)
    if include_length_term:
        score += float(weights.weights["length"]) * math.log(len(normalized))
    return score


def normalize(score: float, weights: CatGranuleWeights) -> float:
    return (score - weights.normalization_mean) / weights.normalization_std


def residue_profile(sequence: str, weights: CatGranuleWeights) -> list[float]:
    """Sliding-window mean over the interior residues.

    A ``residue_window_size``-residue window (radius 25 for the default 51)
    is averaged around each center in the interior span, which yields
    ``len(L) - (residue_window_size - 1)`` values. The first profile value
    maps to sequence position ``offset 24`` (0-indexed), matching the
    PhaSePred web archive residue arrays (FEATURES.md §7). Values are
    un-normalized.
    """
    radius = (weights.residue_window_size - 1) // 2
    if radius <= 0:
        raise CatGranuleInputError("residue_window_size must be at least 3")
    scores = residue_scores(sequence, weights)
    length = len(scores)
    interior_start = radius - 1
    interior_stop = length - radius - 1
    if interior_stop <= interior_start:
        return []
    values: list[float] = []
    for center in range(interior_start, interior_stop):
        lo = max(0, center - radius)
        hi = min(length - 1, center + radius) + 1
        window = scores[lo:hi]
        values.append(sum(window) / len(window))
    return values


def granule_strength(normalized_profile: list[float]) -> float:
    """Equation 3: mean of the positive (Z-normalized) profile values."""
    positive_values = [value for value in normalized_profile if value > 0]
    if not positive_values:
        return 0.0
    return sum(positive_values) / len(positive_values)


def _score_distilled(
    sequence: str,
    weights: CatGranuleWeights,
    *,
    normalized: bool,
) -> dict[str, float | list[float]]:
    """Distilled (surrogate) route scoring.

    - ``single``: the XGBoost surrogate prediction of the web-archive
      catGRANULE value (calendar-absolute caliber).
    - ``residue``: the paper-formula 51-window interior profile, translated so
      its mean equals ``single`` (per-protein affine alignment; preserves the
      well-validated formula geometry while anchoring it to the surrogate).
      With ``normalized=False`` the un-aligned raw formula profile is
      returned instead.
    """
    from catgranule import surrogate

    normalize_sequence(sequence, weights.scales)
    engine = surrogate.load_surrogate(weights)
    single = engine.predict_single(sequence)
    raw_profile = residue_profile(sequence, weights)
    if normalized:
        if raw_profile:
            shift = single - (sum(raw_profile) / len(raw_profile))
            profile: list[float] = [value + shift for value in raw_profile]
        else:
            profile = []
        return {
            "single": single,
            "residue": profile,
            "granule_strength": granule_strength(profile),
        }
    return {
        "single": single,
        "residue": raw_profile,
        "granule_strength": granule_strength(raw_profile),
    }


def score_sequence(
    sequence: str,
    *,
    normalized: bool = True,
    weights: CatGranuleWeights | None = None,
    route: str | None = None,
) -> dict[str, float | list[float]]:
    """Score a single protein sequence.

    Returns ``{"single": float, "residue": list[float], "granule_strength":
    float}``.

    - Formula routes (``paper`` / ``legacy-unaudited``): ``single`` is the
      Z-normalized Equation-2 score by default; ``residue`` is the interior
      sliding-window profile (length ``L-50`` for the default window) with the
      same normalization; ``granule_strength`` is Equation 3.
    - Distilled route (default): ``single`` is the surrogate prediction of the
      PhaSePred web-archive value; ``residue`` is the paper-formula profile
      aligned to ``single`` (mean-equal); ``granule_strength`` is the mean of
      its positive values. With ``normalized=False`` the underlying raw
      paper-formula single and un-aligned profile are returned.
    """
    resolved = weights or load_weights(route=route)
    if resolved.is_distilled:
        return _score_distilled(sequence, resolved, normalized=normalized)
    raw = raw_score(sequence, resolved)
    raw_profile = residue_profile(sequence, resolved)
    if normalized:
        single_value = normalize(raw, resolved)
        norm_profile = [normalize(value, resolved) for value in raw_profile]
    else:
        single_value = raw
        norm_profile = raw_profile
    return {
        "single": single_value,
        "residue": norm_profile,
        "granule_strength": granule_strength(norm_profile),
    }


def score_batch(
    records: Iterable[Mapping[str, str]],
    *,
    normalized: bool = True,
    weights: CatGranuleWeights | None = None,
    route: str | None = None,
) -> Iterator[dict[str, float | list[float] | str]]:
    """Score several sequences iteratively.

    ``records`` is an iterable of ``{"accession": str, "sequence": str}``
    dicts (partitioned batches are handled by the caller). Each yielded
    dict has the ``score_sequence`` shape plus ``accession``. The iteration
    is sequential, single-threaded and deterministic.
    """
    resolved = weights or load_weights(route=route)
    for record in records:
        try:
            accession = str(record["accession"])
            sequence = str(record["sequence"])
        except KeyError as exc:
            raise CatGranuleInputError(
                f"score_batch records must be {{'accession': str, 'sequence': str}}: "
                f"missing {exc.args[0]!r} in {record!r}"
            ) from exc
        scored = score_sequence(sequence, normalized=normalized, weights=resolved)
        yield {"accession": accession, **scored}


def normalize_sequence(sequence: str, scales: Mapping[str, object]) -> str:
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise CatGranuleInputError("Cannot score an empty sequence")
    unknown = sorted({residue for residue in normalized if residue not in scales})
    if unknown:
        raise CatGranuleInputError(
            f"Unsupported residues in sequence: {', '.join(unknown)}"
        )
    return normalized