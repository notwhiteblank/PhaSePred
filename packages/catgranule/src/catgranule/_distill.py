"""Shipped 181-dim feature map for the distilled surrogate (route=distilled).

This module is the package-side port of ``packages/catgranule/scripts/_s4_numeric.py``
(the exact vectorized builder used at distillation training time). It must stay
bit-identical to that script; ``tests/test_distilled.py`` asserts equal outputs
on representative sequences.

The map is pure sequence + numpy (no external tools):

- per-residue (20, 6) scale values from mmc1 Table S4
  (``weights.scales``) and the Bolognesi 2016 published coefficients
  (``weights.weights``);
- centered-window summary statistics (mean/p25/p50/p75/std) at window radii
  3, 7, 13, 25 for all six scales  -> 4 * 5 * 6 = 120 features;
- profile statistics of the Equation-1 residue score at half-window radii
  3, 12, 25  (11 quantile/stat values + fraction positive + mean positive)
  -> 3 * 13 = 39 features;
- amino-acid composition (20 features);
- length and log-length (2 features);
- total 181.
"""

from __future__ import annotations

import math

import numpy as np

from catgranule.weights import SCALE_COLUMNS, CatGranuleWeights

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}
_STANDARD_AA = frozenset(AA_ORDER)

_FEATURE_DIM = 181


def encode_bytes(sequence: str) -> np.ndarray:
    return np.array([_AA_INDEX[ch] for ch in sequence], dtype=np.int32)


def scale_matrix(weights: CatGranuleWeights) -> np.ndarray:
    matrix = np.zeros((20, 6))
    for k, column in enumerate(SCALE_COLUMNS):
        for aa, row in weights.scales.items():
            matrix[_AA_INDEX[aa], k] = row[column]
    return matrix


def coefficients(weights: CatGranuleWeights) -> np.ndarray:
    return np.array([weights.weights[c] for c in SCALE_COLUMNS])


def window_mean(values: np.ndarray, radius: int) -> np.ndarray:
    """(L,6) -> (L,6) truncated-centered window means for a given radius."""
    length, cols = values.shape
    cm = np.zeros((length + 1, cols))
    np.cumsum(values, axis=0, out=cm[1:])
    centers = np.arange(length)
    lo = np.maximum(0, centers - radius)
    hi = np.minimum(length, centers + radius + 1)
    width = (hi - lo).astype(float)[:, None]
    return (cm[hi] - cm[lo]) / width


def residue_scores(
    sequence: str, scales: np.ndarray, params: np.ndarray, radius: int = 3
) -> np.ndarray:
    x = encode_bytes(sequence)
    return window_mean(scales[x], radius) @ params


def window_summary(values: np.ndarray, radius: int) -> list[np.ndarray]:
    wm = window_mean(values, radius)
    return [
        np.mean(wm, axis=0),
        np.quantile(wm, 0.25, axis=0),
        np.quantile(wm, 0.50, axis=0),
        np.quantile(wm, 0.75, axis=0),
        np.std(wm, axis=0),
    ]


def profile_stats(g: np.ndarray, radius: int) -> np.ndarray:
    length = len(g)
    nan = float("nan")
    if length <= 2 * radius:
        return np.full(13, nan)
    cum = np.concatenate([[0.0], np.cumsum(g)])
    centers = np.arange(radius - 1, length - radius - 1)
    lo = np.maximum(0, centers - radius)
    hi = np.minimum(length, centers + radius + 1)
    width = (hi - lo).astype(float)
    prof = (cum[hi] - cum[lo]) / width
    base = [
        np.nanmean(prof),
        np.nanstd(prof),
        np.quantile(prof, 0.05),
        np.quantile(prof, 0.10),
        np.quantile(prof, 0.25),
        np.quantile(prof, 0.50),
        np.quantile(prof, 0.75),
        np.quantile(prof, 0.90),
        np.quantile(prof, 0.95),
        np.nanmin(prof),
        np.nanmax(prof),
    ]
    frac_pos = float(np.mean(prof > 0))
    mean_pos = float(np.nanmean(prof[prof > 0])) if np.any(prof > 0) else 0.0
    return np.array(base + [frac_pos, mean_pos])


def featurize(sequence: str, weights: CatGranuleWeights) -> np.ndarray:
    """(L,) -> (181,) distilled feature vector (NaN regions filled with 0, as
    at training time)."""
    x = encode_bytes(sequence)
    length = x.shape[0]
    scales = scale_matrix(weights)
    params = coefficients(weights)
    vals = scales[x]
    parts: list[np.ndarray] = []
    for radius in (3, 7, 13, 25):
        parts.extend(window_summary(vals, radius))
    g = residue_scores(sequence, scales, params, 3)
    parts.append(profile_stats(g, 3))
    parts.append(profile_stats(g, 12))
    parts.append(profile_stats(g, 25))
    comp = np.zeros(20)
    np.add.at(comp, x, 1)
    comp = comp / length
    parts.append(comp)
    parts.append(np.array([length, math.log(length)], dtype=float))
    vec = np.concatenate(parts)
    assert vec.shape[0] == _FEATURE_DIM
    return np.nan_to_num(vec, nan=0.0, posinf=0.0, neginf=0.0)


def check_match(sequence: str, weights: CatGranuleWeights) -> bool:
    """Sanity helper (used by tests): ensure standard residues only."""
    return set(sequence) <= _STANDARD_AA and len(sequence) > 0


__all__ = [
    "AA_ORDER",
    "_STANDARD_AA",
    "coefficients",
    "encode_bytes",
    "featurize",
    "profile_stats",
    "residue_scores",
    "scale_matrix",
    "window_mean",
    "window_summary",
]