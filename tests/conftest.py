"""Shared fixtures and subprocess helpers for the phasepred test suite.

Centralises three concerns so individual tests do not reinvent them:

- repository / fixture locations (``repo_root``, ``fixtures_dir``);
- tool and private-data availability probes (``tool_available``,
  ``require_tool``, ``private_data_available``) so "missing tool -> SKIP with
  an explicit reason" is structural rather than a per-test convention;
- ``run_cli``, which invokes the installed console script the way a user
  would, pinning ``PHASEPRED_PYTHON`` to the interpreter running the tests
  (Ruling 1) and never the silent-no-op ``python -m phasepred.cli`` form
  (Ruling 5).
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# The six predict-required manifest contracts, in registry order. ``predict``'s
# preflight aborts on the first one that is missing, so a test driving the
# ``predict`` command needs the whole set. catGRANULE is a package import
# (present whenever phasepred is installed) and is not in this list.
PREDICT_REQUIRED_TOOLS: tuple[str, ...] = (
    "seg",
    "plaac",
    "pscore",
    "espritz",
    "deepcoil",
    "localcider",
)

# The two non-contract data rows of the 9-row check-tools registry. They have no
# tools/<KEY>/manifest.toml (DeepPhase ships bundled; PhosphoSitePlus is
# registration-gated), so ``require_tool`` cannot see them and the probes below
# resolve them the way their consumers do.
_OPTIONAL_DATA_PROBES: dict[str, str] = {
    "deepphase": "deepphase_table_path",
    "phosphositeplus": "phosphosite_path",
}

# The tests escalated during E8's preparation (D43 finding 1) plus the one
# further case that only a genuinely isolated public-checkout run exposed: each
# needs an external tool or private data that is gitignored, so each carries a
# conditional skip guard. ``test_ci_workflow`` proves they SKIP (not fail) in a
# tool-less tree, and ``test_optional_tool_skips`` proves they still PASS in a
# provisioned one, so the guard cannot become permanent silence.
ESCALATION_TOOL_TESTS: tuple[str, ...] = (
    "tests/test_cli.py::test_features_recomputed_requires_external_deepcoil",
    "tests/test_cli.py::test_features_recomputed_csv_source_requires_catgranule_features",
    "tests/test_cli.py::test_predict_no_product_flag_uses_single_product",
    "tests/test_cli.py::test_predict_stale_warning_suppressed_fresh_manifest",
    "tests/test_tool_features.py::test_lookup_degrades_to_nan",
    "tests/test_tool_packages.py::test_contract_run_check_all_six",
    "tests/test_tool_packages.py::test_contract_install_check_all_six",
    "tests/test_tool_packages.py::test_contract_run_produces_stdout[pscore]",
    "tests/test_tool_packages.py::test_contract_run_produces_stdout[espritz]",
    "tests/test_tool_packages.py::test_check_tools_cli_strict_exits_zero",
)

# Gitignored private inputs (E5 golden fixtures are derived from these). The
# name -> repo-relative path table is explicit so a typo fails loudly instead
# of silently reporting "absent".
_PRIVATE_DATA: dict[str, tuple[str, ...]] = {
    "human_reviewed": (
        "data",
        "raw",
        "external",
        "phasepred_web",
        "human_reviewed.json",
    ),
    "human_compact": (
        "data",
        "raw",
        "external",
        "phasepred_web",
        "human_compact.jsonl",
    ),
    "uniprot_cache": ("data", "interim", "uniprot_cache.jsonl"),
}


@pytest.fixture
def repo_root() -> Path:
    """The repository root (the checkout or the wheel's installed root)."""
    return REPO_ROOT


@pytest.fixture
def fixtures_dir() -> Path:
    """The committed test-fixture directory (``tests/fixtures``)."""
    return REPO_ROOT / "tests" / "fixtures"


def _resolve_runner() -> list[str]:
    """Return the argv prefix that runs the phasepred CLI.

    Ruling 5: ``python -m phasepred.cli`` is a silent no-op (exit 0, empty
    stdout, no file written), so it is never the default. Preference:

    1. ``$PHASEPRED_RUNNER`` — the tests parameterise over this; split with
       :func:`shlex.split` so it may carry arguments.
    2. the ``phasepred`` console script on ``PATH``.
    3. the console script beside the running interpreter (``sys.executable``),
       which is where a venv/pixi environment installs it.

    The chosen form's output is pinned by
    ``test_cwd.py::test_cli_runner_produces_parseable_output``.
    """
    override = os.environ.get("PHASEPRED_RUNNER")
    if override:
        return shlex.split(override)
    found = shutil.which("phasepred")
    if found:
        return [found]
    sibling = Path(sys.executable).with_name("phasepred")
    if sibling.is_file():
        return [str(sibling)]
    raise RuntimeError(
        "no phasepred runner found: set PHASEPRED_RUNNER or install the "
        "console script next to sys.executable"
    )


def run_cli(
    *args: str,
    env: dict[str, str] | None = None,
    cwd: str | Path | None = None,
    timeout: float = 300,
) -> subprocess.CompletedProcess[str]:
    """Run the phasepred CLI as a subprocess and return the completed process.

    ``PHASEPRED_PYTHON`` is injected (Ruling 1) so contract ``run`` scripts
    resolve the same interpreter the parent runs under; ``env`` merges on top.
    The invocation form comes from :func:`_resolve_runner` (Ruling 5).
    """
    child_env = {**os.environ, "PHASEPRED_PYTHON": sys.executable}
    if env:
        child_env.update(env)
    return subprocess.run(
        [*_resolve_runner(), *args],
        capture_output=True,
        text=True,
        env=child_env,
        cwd=cwd,
        timeout=timeout,
    )


def tool_available(key: str) -> bool:
    """Whether tool ``key`` resolves *and* passes its ``run --check`` probe.

    ``key`` is the registry key (e.g. ``"espritz"``); the probe mirrors
    :func:`phasepred.tool_paths.check_all_tools` for that one row.
    """
    try:
        from phasepred.tool_paths import _probe_usable, resolve_status

        status = resolve_status(key)
    except KeyError:
        return False
    if not status.ok or status.path is None:
        return False
    if status.source == "path":
        return True
    return _probe_usable(status.path)


def require_tool(key: str) -> None:
    """``pytest.skip`` with an explicit reason unless tool ``key`` is usable."""
    if tool_available(key):
        return
    reason = f"no manifest registered for tool '{key}'"
    try:
        from phasepred.tool_paths import resolve_status

        status = resolve_status(key)
        if status.ok and status.path is not None:
            # tool_available() only rejects this row because the contract's
            # ``run --check`` probe failed, so say that instead of repeating
            # the install hint (which would be misleading here).
            reason = f"{status.path!r} failed its --check probe"
        else:
            reason = status.install_hint or f"{key} did not resolve"
    except KeyError:
        pass
    pytest.skip(f"tool '{key}' unavailable: {reason}")


def require_predict_requirements() -> None:
    """Skip unless every predict-required contract resolves and probes usable.

    ``predict`` aborts on the first missing required tool, so any test that
    drives that command must see the full set. Composes :func:`require_tool` so
    the skip reason stays the single per-tool reason that helper produces.
    """
    for key in PREDICT_REQUIRED_TOOLS:
        require_tool(key)


def optional_data_available(name: str) -> bool:
    """Whether a non-contract data row (``deepphase``/``phosphositeplus``) exists.

    Resolved through the same lazy accessor its consumer uses, so setting
    ``PHASEPRED_DATA_ROOT`` moves this probe exactly as it moves the feature.
    """
    from phasepred import tools

    try:
        accessor = _OPTIONAL_DATA_PROBES[name]
    except KeyError as e:
        raise KeyError(
            f"unknown optional data {name!r}; known: {sorted(_OPTIONAL_DATA_PROBES)}"
        ) from e
    return bool(getattr(tools, accessor)().is_file())


def require_optional_data(name: str) -> None:
    """``pytest.skip`` with an explicit reason unless the data row is installed."""
    if optional_data_available(name):
        return
    pytest.skip(
        f"optional data '{name}' is not installed in this checkout "
        f"(run the matching tools/<TOOL>/install.sh or set PHASEPRED_DATA_ROOT)"
    )


def private_data_available(name: str) -> bool:
    """Whether a gitignored private input ``name`` is present.

    ``name`` is one of :data:`_PRIVATE_DATA` (``human_reviewed``,
    ``human_compact``, ``uniprot_cache``). Resolved under the active data
    root, so it honours ``PHASEPRED_DATA_ROOT``.
    """
    try:
        parts = _PRIVATE_DATA[name]
    except KeyError as e:
        raise KeyError(
            f"unknown private data {name!r}; known: {sorted(_PRIVATE_DATA)}"
        ) from e
    from phasepred.data import data_path

    return data_path(*parts).is_file()


def installed_psp_version_marker() -> str | None:
    """Return the locally installed PhosphoSitePlus version marker.

    Reads the first line of the gunzipped dataset (Ruling 2). The
    ``PHASEPRED_PSP_VERSION_MARKER`` override exists so the mismatch->SKIP
    branch can be exercised without touching the real data; it is the reverse
    proof for gate T3-G4 and is never set in normal runs.
    """
    override = os.environ.get("PHASEPRED_PSP_VERSION_MARKER")
    if override:
        return override
    from phasepred.tools import phosphosite_version_marker

    return phosphosite_version_marker()


def phos_freq_matches(
    frozen: float,
    actual: float,
    *,
    frozen_marker: str,
    installed_marker: str | None = None,
) -> None:
    """Compare a local ``Phos freq`` value to its frozen gold (Ruling 2).

    Marker match -> exact equality; marker mismatch -> ``pytest.skip`` with a
    stated reason (never FAIL, never PASS). Task 4's functional matrix and
    E6/E9's smoke harness call this so a user whose PhosphoSitePlus has rolled
    forward gets a SKIP, not a spurious failure.
    """
    marker = installed_marker if installed_marker is not None else installed_psp_version_marker()
    if marker != frozen_marker:
        pytest.skip(
            f"PhosphoSitePlus version mismatch: local marker {marker!r} != "
            f"frozen {frozen_marker!r}; re-freeze the Phos freq gold after a "
            "PSP update (Ruling 2)"
        )
    assert actual == frozen, f"Phos freq {actual!r} != frozen {frozen!r}"


@pytest.fixture
def psp_env() -> SimpleNamespace:
    """Expose the Ruling-2 Phos freq gate and the local PSP marker probe."""
    return SimpleNamespace(
        installed_marker=installed_psp_version_marker,
        matches=phos_freq_matches,
    )
