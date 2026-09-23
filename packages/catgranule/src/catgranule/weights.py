"""Weight-artifact loading and provenance validation.

The scoring parameters are externalized as replaceable JSON artifacts inside
the package so they can be swapped as a unit. Two classes of artifacts exist
(selection via ``CATGRANULE_WEIGHTS`` env var, the ``route`` argument, or the
``--route`` CLI flag):

- ``formula`` artifacts (``mode="formula"``): the Bolognesi 2016 closed-form
  recipe — Equation 1-2 coefficients, the mmc1 Table S4 residue scales and
  the yeast-proteome Z-normalization constants. Routes ``paper`` (S4
  recomputed Z constants) and ``legacy-unaudited`` (pre-S4 reconstruction)
  are formula artifacts.
- ``distilled`` artifact (``mode="distilled"``): the S4-D21 XGBoost surrogate
  that predicts the PhaSePred web-archive catGRANULE ``single`` directly. It
  embeds its own copy of the formula scales/coefficients (its feature basis)
  plus the path of the serialized XGBoost model.

Each artifact carries a ``provenance`` block describing how it was produced
and whether it has passed any gate.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:  # importlib.resources.files is available on Python 3.9+
    from importlib.resources import files as _resources_files
except ImportError:  # pragma: no cover - Python <3.9 fallback
    _resources_files = None  # type: ignore[assignment]

DEFAULT_ARTIFACT = "catgranule_distilled_weights.json"
ENV_OVERRIDE = "CATGRANULE_WEIGHTS"

VALID_ROUTES = frozenset({"paper", "distilled", "legacy-unaudited"})
SCALE_COLUMNS = ("rc", "rn", "dc", "dn", "pfg", "prg")

ROUTE_ARTIFACTS: dict[str, str] = {
    "distilled": "catgranule_distilled_weights.json",
    "paper": "catgranule_paper_weights.json",
    "legacy-unaudited": "catgranule_v1_weights.json",
}


class CatGranuleWeightsError(RuntimeError):
    """Raised when a weights artifact is missing or invalid."""


@dataclass(frozen=True)
class Provenance:
    route: str
    source: str
    date: str
    gate: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "route": self.route,
            "source": self.source,
            "date": self.date,
            "gate": self.gate,
        }


@dataclass(frozen=True)
class CatGranuleWeights:
    version: str
    weights: dict[str, float]
    scales: dict[str, dict[str, float]]
    normalization_mean: float
    normalization_std: float
    normalization_n_proteins: int
    window_radius: int
    residue_window_size: int
    provenance: Provenance
    mode: str = "formula"
    model_path: Path | None = None
    feature_spec: dict[str, Any] | None = None

    @property
    def is_legacy_unaudited(self) -> bool:
        return self.provenance.route == "legacy-unaudited"

    @property
    def is_distilled(self) -> bool:
        return self.mode == "distilled"

    def to_json(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "weights": self.weights,
            "scales": self.scales,
            "normalization": {
                "mean": self.normalization_mean,
                "std": self.normalization_std,
                "n_proteins": self.normalization_n_proteins,
            },
            "window_radius": self.window_radius,
            "residue_window_size": self.residue_window_size,
            "mode": self.mode,
            "model_file": self.model_path.name if self.model_path else None,
            "feature_spec": self.feature_spec,
            "provenance": self.provenance.to_json(),
        }


def _packaged_artifact_path(name: str) -> Path:
    if _resources_files is not None:
        resource = _resources_files("catgranule") / "weights" / name
        return Path(resource)  # type: ignore[arg-type]
    return Path(__file__).with_name("weights") / name


def _default_artifact_path() -> Path:
    return _packaged_artifact_path(DEFAULT_ARTIFACT)


def resolve_artifact_path(route: str | None = None, path: Path | None = None) -> Path:
    """Resolve the artifact path: explicit ``path`` > env override >
    packaged ``route`` > default."""
    if path is not None:
        return path.expanduser().resolve()
    override = os.environ.get(ENV_OVERRIDE)
    if override:
        return Path(override).expanduser().resolve()
    if route is not None:
        if route not in ROUTE_ARTIFACTS:
            raise CatGranuleWeightsError(
                f"Unknown weights route {route!r}; expected one of {sorted(ROUTE_ARTIFACTS)}"
            )
        return _packaged_artifact_path(ROUTE_ARTIFACTS[route])
    return _default_artifact_path()


def _read_artifact(path: Path) -> dict[str, Any]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise CatGranuleWeightsError(
            f"Weights artifact not found: {path}. "
            f"Set {ENV_OVERRIDE} to a valid artifact or reinstall the package."
        ) from exc
    if not isinstance(raw, dict):
        raise CatGranuleWeightsError(f"Weights artifact {path} is not a JSON object")
    return raw


def _load_raw() -> dict[str, Any]:
    """Back-compat helper: read the effective artifact raw dict."""
    return _read_artifact(resolve_artifact_path())


def _parse_provenance(raw: dict[str, Any]) -> Provenance:
    provenance = raw.get("provenance")
    if not isinstance(provenance, dict):
        raise CatGranuleWeightsError("Weights artifact missing 'provenance' block")
    route = provenance.get("route")
    if route not in VALID_ROUTES:
        raise CatGranuleWeightsError(
            f"Unknown provenance route {route!r}; expected one of {sorted(VALID_ROUTES)}"
        )
    source = provenance.get("source")
    date = provenance.get("date")
    if not isinstance(source, str) or not isinstance(date, str):
        raise CatGranuleWeightsError("provenance 'source' and 'date' must be strings")
    gate = provenance.get("gate")
    return Provenance(route=route, source=source, date=date, gate=gate)


def load_weights(
    route: str | None = None, path: Path | None = None
) -> CatGranuleWeights:
    """Load a weights artifact (default distilled; ``route``/``path``/env override)."""
    resolved_path = resolve_artifact_path(route, path)
    raw = _read_artifact(resolved_path)

    weights = raw.get("weights")
    scales = raw.get("scales")
    if not isinstance(weights, dict) or not isinstance(scales, dict):
        raise CatGranuleWeightsError("Weights artifact missing 'weights' or 'scales'")

    mode = str(raw.get("mode", "formula"))
    if mode not in ("formula", "distilled"):
        raise CatGranuleWeightsError(f"Unknown artifact mode {mode!r}")

    normalization = raw.get("normalization")
    model_path: Path | None = None
    if mode == "distilled":
        model_file = raw.get("model_file")
        if not isinstance(model_file, str) or not model_file:
            raise CatGranuleWeightsError(
                "Distilled artifact requires a 'model_file' entry"
            )
        model_path = resolved_path.parent / model_file
        if not model_path.exists():
            raise CatGranuleWeightsError(
                f"Distilled model file not found next to the artifact: {model_path}. "
                f"Set {ENV_OVERRIDE}/route to a complete artifact directory or reinstall."
            )
        if not isinstance(normalization, dict) or "mean" not in normalization:
            raise CatGranuleWeightsError(
                "Distilled artifact still carries a normalization block (paper-formula "
                "basis); missing or incomplete."
            )
    elif not isinstance(normalization, dict) or (
        "mean" not in normalization or "std" not in normalization
    ):
        raise CatGranuleWeightsError("Weights artifact missing normalization mean/std")

    unknown_aas = set(scales) - set("ACDEFGHIKLMNPQRSTVWY")
    if unknown_aas:
        raise CatGranuleWeightsError(
            f"Weights artifact contains unknown residue entries: {sorted(unknown_aas)}"
        )
    for aa, values in scales.items():
        missing = set(SCALE_COLUMNS) - set(values)
        if missing:
            raise CatGranuleWeightsError(
                f"Scale for residue {aa!r} is missing columns: {sorted(missing)}"
            )

    feature_spec = raw.get("feature_spec")
    if feature_spec is not None and not isinstance(feature_spec, dict):
        raise CatGranuleWeightsError("feature_spec must be a JSON object")

    return CatGranuleWeights(
        version=str(raw.get("version", "")),
        weights={col: float(weights[col]) for col in SCALE_COLUMNS + ("length",)},
        scales={
            aa: {col: float(values[col]) for col in SCALE_COLUMNS}
            for aa, values in scales.items()
        },
        normalization_mean=float(normalization["mean"]) if "mean" in normalization else 0.0,
        normalization_std=float(normalization["std"]) if "std" in normalization else 1.0,
        normalization_n_proteins=int(normalization.get("n_proteins", 0)),
        window_radius=int(raw.get("window_radius", 3)),
        residue_window_size=int(raw.get("residue_window_size", 51)),
        provenance=_parse_provenance(raw),
        mode=mode,
        model_path=model_path,
        feature_spec=feature_spec,
    )