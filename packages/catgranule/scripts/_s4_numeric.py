"""Shared vectorized feature builders for S4 catGRANULE retraining.

These mirror ``catgranule.scoring`` exactly (verified bit-identical to the
package implementation for raw score and residue profile, diff < 2e-14) but
use numpy so that the several-thousand-protein evaluation loops of the S4
gate run in seconds instead of minutes. All numerics are deterministic and
written to match the package computation of Equation 1-2 (truncated centered
heptapeptide windows, mean over the protein, `+ 0.25*log(L)`).
"""

from __future__ import annotations

import math

import numpy as np

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}
_STANDARD_AA = frozenset(AA_ORDER)


def encode_bytes(sequence: str) -> np.ndarray:
    return np.array([_AA_INDEX[ch] for ch in sequence], dtype=np.int32)


def make_scale_matrix(weights) -> np.ndarray:
    """Return the (20, 6) residue-scale matrix from a CatGranuleWeights object."""
    from catgranule.weights import SCALE_COLUMNS

    matrix = np.zeros((20, 6))
    for k, column in enumerate(SCALE_COLUMNS):
        for aa, row in weights.scales.items():
            matrix[_AA_INDEX[aa], k] = row[column]
    return matrix


def paper_coefficients(weights) -> np.ndarray:
    from catgranule.weights import SCALE_COLUMNS

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
    sequence: str, scale_matrix: np.ndarray, coefficients: np.ndarray, radius: int = 3
) -> np.ndarray:
    x = encode_bytes(sequence)
    return window_mean(scale_matrix[x], radius) @ coefficients


def raw_score(
    sequence: str, scale_matrix: np.ndarray, coefficients: np.ndarray,
    a_length: float = 0.25, radius: int = 3,
) -> float:
    length = len(sequence)
    g = residue_scores(sequence, scale_matrix, coefficients, radius)
    return float(g.mean() + a_length * math.log(length))


def profile(
    sequence: str, scale_matrix: np.ndarray, coefficients: np.ndarray, radius: int = 25
) -> np.ndarray:
    """Interior 51-window mean profile (same geometry as ``residue_profile``)."""
    g = residue_scores(sequence, scale_matrix, coefficients, 3)
    length = len(g)
    if length <= 2 * radius:
        return np.empty(0)
    cum = np.concatenate([[0.0], np.cumsum(g)])
    centers = np.arange(radius - 1, length - radius - 1)
    lo = np.maximum(0, centers - radius)
    hi = np.minimum(length, centers + radius + 1)
    width = (hi - lo).astype(float)
    return (cum[hi] - cum[lo]) / width


def window_summary(values: np.ndarray, radius: int) -> list[np.ndarray]:
    """Summary of radius-banded window means: mean/p25/p50/p75/std per scale."""
    wm = window_mean(values, radius)
    return [
        np.mean(wm, axis=0),
        np.quantile(wm, 0.25, axis=0),
        np.quantile(wm, 0.50, axis=0),
        np.quantile(wm, 0.75, axis=0),
        np.std(wm, axis=0),
    ]


def profile_stats(g: np.ndarray, radius: int) -> np.ndarray:
    """Summary stats of a g-based profile at a given half-window radius."""
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


def exposure_features(sequence: str, radius: int = 3) -> np.ndarray:
    """Per-residue effective-exposure weights folded into 20 AA features.

    For residue ``i`` the contribution of every residue inside its centered
    window is ``1/width``; summing over positions yields, per amino-acid, an
    exposure weight. Combined with ``log(L)`` this spans the exact linear
    functional form of Equation 1-2 with free effective coefficients.
    """
    length = len(sequence)
    w = np.zeros(length)
    for i in range(length):
        lo = max(0, i - radius)
        hi = min(length, i + radius + 1)
        w[lo:hi] += 1.0 / (hi - lo)
    feat = np.zeros(20)
    x = encode_bytes(sequence)
    np.add.at(feat, x, w)
    feat /= length
    return feat


def featurize_distill(
    sequence: str, scale_matrix: np.ndarray, coefficients: np.ndarray
) -> np.ndarray:
    """Rich (sequence -> vector) feature map used for the distilled model."""
    x = encode_bytes(sequence)
    length = x.shape[0]
    vals = scale_matrix[x]
    parts: list[np.ndarray] = []
    for radius in (3, 7, 13, 25):
        parts.extend(window_summary(vals, radius))
    g = residue_scores(sequence, scale_matrix, coefficients, 3)
    parts.append(profile_stats(g, 3))
    parts.append(profile_stats(g, 12))
    parts.append(profile_stats(g, 25))
    comp = np.zeros(20)
    np.add.at(comp, x, 1)
    comp = comp / length
    parts.append(comp)
    parts.append(np.array([length, math.log(length)], dtype=float))
    return np.concatenate(parts)