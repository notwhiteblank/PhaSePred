"""catGRANULE phase-separation propensity scoring.

Public API (used by ``phasepred``):

- ``score_sequence(sequence, *, normalized=True, route=None)`` -> dict with
  keys ``single`` (protein-level score), ``residue`` (interior
  sliding-window profile, length ``L-50``), ``granule_strength``.
- ``score_batch(records, *, normalized=True, route=None)`` -> iterator over
  the same shape, where ``records`` is an iterable of ``{"accession": str,
  "sequence": str}`` dicts; each result adds an ``accession`` key.

Weights are externalized as replaceable artifacts (see
``catgranule.weights``); ``route`` selects among the packaged artifacts
(``distilled`` default, ``paper``, ``legacy-unaudited``) and
``$CATGRANULE_WEIGHTS`` or ``--weights`` overrides with a file path:

- ``distilled`` (default): XGBoost surrogate that predicts the PhaSePred
  web-archive catGRANULE ``single`` (pure sequence features, no external
  tools).
- ``paper``: Bolognesi 2016 closed-form recipe with S4-recomputed yeast Z
  constants (alternate artifact).
- ``legacy-unaudited``: the pre-S4 reconstruction (formula compatibility).
"""

from __future__ import annotations

from catgranule.scoring import (
    CatGranuleInputError,
    score_batch,
    score_sequence,
)
from catgranule.weights import load_weights

__version__ = "1.0.1"
__all__ = [
    "CatGranuleInputError",
    "load_weights",
    "score_batch",
    "score_sequence",
    "__version__",
]