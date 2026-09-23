"""Canonical root / data-file resolution for the phasepred package.

Every repo-local or shipped data path used by the package goes through the
helpers here so that the CLI and API work from any working directory:

- ``PHASEPRED_DATA_ROOT`` overrides where repo-local data lives
  (``data/``, ``tools/``, ``products/``, ``.external_envs/``). When the
  package is pip-installed without a checkout, point this at a directory
  laid out like the repository root.
- Without an override, the repository root is inferred from this module's
  location (package is installed from ``src/phasepred/`` inside a checkout,
  so ``parents[2]`` is the repo root).
- Package-bundled resources (JSON calibration constants, etc.) are resolved
  relative to the module itself.
"""

from __future__ import annotations

import os
from pathlib import Path


def repo_root() -> Path:
    """The repository root, inferred from the module location."""
    return Path(__file__).resolve().parents[2]


def data_root() -> Path:
    """The directory that holds repo-local data (override-able).

    Preference: ``PHASEPRED_DATA_ROOT`` env var > inferred repository root.
    """
    override = os.environ.get("PHASEPRED_DATA_ROOT")
    if override:
        return Path(override).expanduser().resolve()
    return repo_root()


def data_path(*parts: str) -> Path:
    """Resolve a path under the data root (e.g. ``data_path("data", "raw", ...)``)."""
    return data_root().joinpath(*parts)


def bundled_path(*parts: str) -> Path:
    """Resolve a package-bundled resource (inside ``src/phasepred/``)."""
    return Path(__file__).parent.joinpath(*parts)


def models_root() -> Path:
    """Default Product A artifacts root.

    Preference: ``PHASEPRED_DATA_ROOT`` override (external data root laid out
    like the repo) > package-bundled models (``src/phasepred/data/models``,
    shipped in the wheel) > inferred repository ``products/`` layout (dev
    checkout fallback).
    """
    override = os.environ.get("PHASEPRED_DATA_ROOT")
    if override:
        return Path(override).expanduser().resolve().joinpath(
            "products", "A_paper_split_recomputed", "models"
        )
    bundled = bundled_path("data", "models")
    if bundled.is_dir():
        return bundled
    return data_path("products", "A_paper_split_recomputed", "models")


# Optional external data is laid out flat under a user-writable data root
# (``<root>/<component>/...``); historical checkouts kept the same components
# under ``data/raw/external/<legacy-name>/...``. The mapping is an explicit
# table rather than string guessing so a renamed component can never silently
# resolve to the wrong legacy directory.
_COMPONENT_LEGACY_NAMES: dict[str, str] = {
    "deepphase": "deepphase",
    "phosphositeplus": "phosphositeplus",
    "pscore": "pscore",
    "espritz": "espritz",
}


def _xdg_data_root() -> Path:
    """``$XDG_DATA_HOME`` (or ``~/.local/share`` when unset), resolved."""
    xdg = os.environ.get("XDG_DATA_HOME")
    if xdg:
        return Path(xdg).expanduser().resolve()
    return Path("~/.local/share").expanduser().resolve()


def user_data_root() -> Path:
    """Preferred user-writable root for optional external data.

    Preference: ``PHASEPRED_DATA_ROOT`` > ``$XDG_DATA_HOME/phasepred`` >
    ``~/.local/share/phasepred``.

    Evaluated on every call (no import-time caching) so a session that sets an
    environment variable after importing :mod:`phasepred.tools` still resolves
    the new location. This function never creates directories.
    """
    override = os.environ.get("PHASEPRED_DATA_ROOT")
    if override:
        return Path(override).expanduser().resolve()
    return _xdg_data_root().joinpath("phasepred")


def resolve_optional_data(component: str, *parts: str) -> tuple[Path, str]:
    """Resolve an optional external-data file and report which tier supplied it.

    Same three-tier probe as :func:`optional_data_path`, but returns
    ``(path, tier)``. ``tier`` is exactly one of:

    - ``"override"`` — a hit under ``$PHASEPRED_DATA_ROOT`` (flat or legacy
      ``data/raw/external/<legacy-name>/`` layout).
    - ``"xdg"`` — a hit under an explicit ``$XDG_DATA_HOME/phasepred``.
    - ``"home"`` — a hit under the ``~/.local/share/phasepred`` fallback
      (only probed when ``XDG_DATA_HOME`` is unset).
    - ``"repo"`` — a hit under ``repo_root()/data/raw/external/<legacy-name>/``
      (the dev-checkout fallback).
    - ``"none"`` — nothing exists; the returned path is the preferred target
      under :func:`user_data_root` (Ruling 7), **not** an existing file.

    Per Ruling 9 this is the tier-reporting counterpart to
    :func:`optional_data_path`; callers that expose provenance to users must
    map ``override``/``xdg``/``home`` to a user-controlled root and ``repo`` to
    the checkout, rather than labelling every hit the same. Never creates
    directories, and evaluates the environment on every call.
    """
    legacy = _COMPONENT_LEGACY_NAMES[component]
    override = os.environ.get("PHASEPRED_DATA_ROOT")
    if override:
        root = Path(override).expanduser().resolve()
        flat = root.joinpath(component, *parts)
        if flat.exists():
            return flat, "override"
        legacy_override = root.joinpath("data", "raw", "external", legacy, *parts)
        if legacy_override.exists():
            return legacy_override, "override"
    xdg_env = os.environ.get("XDG_DATA_HOME")
    if xdg_env:
        xdg_candidate = (
            Path(xdg_env).expanduser().resolve().joinpath("phasepred", component, *parts)
        )
        if xdg_candidate.exists():
            return xdg_candidate, "xdg"
    else:
        home_candidate = (
            Path("~/.local/share").expanduser().resolve().joinpath(
                "phasepred", component, *parts
            )
        )
        if home_candidate.exists():
            return home_candidate, "home"
    repo_candidate = repo_root().joinpath("data", "raw", "external", legacy, *parts)
    if repo_candidate.exists():
        return repo_candidate, "repo"
    return user_data_root().joinpath(component, *parts), "none"


def optional_data_path(component: str, *parts: str) -> Path:
    """Resolve an optional external-data file in three tiers.

    Probed in order for ``<component>`` / ``<parts>``:

    1. ``$PHASEPRED_DATA_ROOT/<component>/<parts>``, also accepting the legacy
       ``$PHASEPRED_DATA_ROOT/data/raw/external/<legacy-name>/<parts>`` layout.
    2. ``$XDG_DATA_HOME/phasepred/<component>/<parts>`` (or
       ``~/.local/share/phasepred/...`` when ``XDG_DATA_HOME`` is unset).
    3. ``repo_root()/data/raw/external/<legacy-name>/<parts>`` — the dev
       checkout fallback that keeps existing local data working.

    Returns the first existing candidate. Per Ruling 7, when none exists it
    returns the target path under :func:`user_data_root` (the preferred write
    root for installers) rather than ``None``. **Callers must therefore call
    ``.exists()`` on the result** — a returned path does not guarantee the data
    is present. This function never creates directories.

    This is a thin wrapper over :func:`resolve_optional_data` (which also
    reports the tier); its signature and behaviour are unchanged.
    """
    return resolve_optional_data(component, *parts)[0]
