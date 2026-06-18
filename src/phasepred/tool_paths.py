"""External tool path resolution for PhaSePred.

Every feature-tool path used by `predict` and the legacy feature pipeline
goes through this module. The resolution cascade is intentionally short:

  1. environment variable  (e.g. PHASEPRED_SEG_BIN, PHASEPRED_ESPRITZ_DIR)
  2. vendored under tools/per-tool/<TOOL>/<entrypoint>
  3. shutil.which on PATH
  4. raise PhaSePredToolNotFound with the install hint string

`check_all_tools()` returns a `ToolStatus` list suitable for the
`phasepred check-tools` CLI command and for early-abort gating inside
`predict`.

The resolver does not invoke any tool. It only reports presence /
absence of files and binaries so callers can decide to proceed or to
print a useful error.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parents[2]
VENDORED_ROOT = REPO_ROOT / "tools" / "per-tool"
WRAPPERS_ROOT = REPO_ROOT / "tools" / "wrappers"


class PhaSePredToolNotFound(RuntimeError):
    """Raised when a required tool cannot be located by the resolver."""


@dataclass
class ToolStatus:
    """Resolution outcome for a single tool."""

    name: str
    ok: bool
    source: str  # "env:<VAR>", "vendored", "path", "built-in", or "missing"
    path: Path | None
    required_for_predict: bool
    install_hint: str = ""


@dataclass
class _ToolSpec:
    """How to look up one tool."""

    name: str
    env_var: str
    vendored_path: Path | None
    path_executables: tuple[str, ...] = ()
    required_for_predict: bool = True
    install_hint: str = ""
    # Optional verifier callable: (path) -> bool
    verifier: Callable[[Path], bool] | None = None
    # Some tools resolve to a directory rather than a single file
    expect_dir: bool = False


def _exists_file(p: Path) -> bool:
    return p.is_file()


def _exists_dir(p: Path) -> bool:
    return p.is_dir()


def _resolve(spec: _ToolSpec) -> ToolStatus:
    """Run the 4-tier cascade for one tool."""

    exists = _exists_dir if spec.expect_dir else _exists_file
    verifier = spec.verifier or exists

    # 1. environment variable
    env_val = os.environ.get(spec.env_var)
    if env_val:
        cand = Path(env_val).expanduser().resolve()
        if verifier(cand):
            return ToolStatus(
                name=spec.name,
                ok=True,
                source=f"env:{spec.env_var}",
                path=cand,
                required_for_predict=spec.required_for_predict,
                install_hint=spec.install_hint,
            )

    # 2. vendored
    if spec.vendored_path is not None and verifier(spec.vendored_path):
        return ToolStatus(
            name=spec.name,
            ok=True,
            source="vendored",
            path=spec.vendored_path,
            required_for_predict=spec.required_for_predict,
            install_hint=spec.install_hint,
        )

    # 3. PATH
    for exe in spec.path_executables:
        hit = shutil.which(exe)
        if hit:
            return ToolStatus(
                name=spec.name,
                ok=True,
                source="path",
                path=Path(hit),
                required_for_predict=spec.required_for_predict,
                install_hint=spec.install_hint,
            )

    # 4. missing
    return ToolStatus(
        name=spec.name,
        ok=False,
        source="missing",
        path=None,
        required_for_predict=spec.required_for_predict,
        install_hint=spec.install_hint,
    )


# ---------------------------------------------------------------------------
# Per-tool specs
# ---------------------------------------------------------------------------

_SEG = _ToolSpec(
    name="SEG",
    env_var="PHASEPRED_SEG_BIN",
    vendored_path=VENDORED_ROOT / "SEG" / "seg",
    path_executables=("seg",),
    required_for_predict=True,
    install_hint=(
        "Build the vendored copy: (cd tools/per-tool/SEG && make), "
        "or set PHASEPRED_SEG_BIN to your seg binary."
    ),
)

_PSCORE = _ToolSpec(
    name="PScore",
    env_var="PHASEPRED_PSCORE_DIR",
    vendored_path=VENDORED_ROOT / "PScore" / "SourceCodeS2",
    expect_dir=True,
    required_for_predict=True,
    install_hint=(
        "PScore source ships vendored under tools/per-tool/PScore/SourceCodeS2/. "
        "If you moved it, set PHASEPRED_PSCORE_DIR to the directory containing "
        "elife_phase_separation_predictor.py."
    ),
    verifier=lambda p: (p / "elife_phase_separation_predictor.py").is_file(),
)

_PLAAC_WRAPPER = WRAPPERS_ROOT / "run_plaac.sh"
_PLAAC = _ToolSpec(
    name="PLAAC",
    env_var="PHASEPRED_PLAAC_JAR",
    vendored_path=VENDORED_ROOT / "PLAAC" / "plaac-master" / "web" / "bin" / "plaac.jar",
    required_for_predict=True,
    install_hint=(
        "PLAAC jar ships vendored at tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar. "
        "Set PHASEPRED_PLAAC_JAR to override. PLAAC also requires a Java runtime on PATH."
    ),
)

_IUPRED3 = _ToolSpec(
    name="IUPred3",
    env_var="PHASEPRED_IUPRED3_DIR",
    vendored_path=VENDORED_ROOT / "IUPred3" / "iupred3",
    expect_dir=True,
    required_for_predict=False,  # optional alternative to ESpritz
    install_hint=(
        "IUPred3 cannot be redistributed (academic license). Accept terms at "
        "https://iupred3.elte.hu/, extract iupred3.tar.gz under "
        "tools/per-tool/IUPred3/iupred3/, or set PHASEPRED_IUPRED3_DIR."
    ),
    verifier=lambda p: (p / "iupred3_lib.py").is_file(),
)

_ESPRITZ = _ToolSpec(
    name="ESpritz",
    env_var="PHASEPRED_ESPRITZ_DIR",
    vendored_path=VENDORED_ROOT / "ESpritz" / "espritz",
    expect_dir=True,
    required_for_predict=True,  # paper-defined IDR feature
    install_hint=(
        "ESpritz cannot be redistributed (academic license). Place espritz.zip "
        "under tools/per-tool/ESpritz/ and run: "
        "(cd tools/per-tool/ESpritz && unzip espritz.zip), "
        "or set PHASEPRED_ESPRITZ_DIR to point to an existing install. "
        "Also requires Perl on PATH."
    ),
    verifier=lambda p: (p / "espritz.pl").is_file(),
)

_DEEPCOIL = _ToolSpec(
    name="DeepCoil",
    env_var="PHASEPRED_DEEPCOIL_ENV",
    vendored_path=REPO_ROOT / ".external_envs" / "deepcoil",
    expect_dir=True,
    required_for_predict=True,
    install_hint=(
        "DeepCoil needs an isolated Python 3.8 conda env. Create it with: "
        "bash tools/install/install_deepcoil_env.sh, "
        "or set PHASEPRED_DEEPCOIL_ENV to an existing prefix."
    ),
    verifier=lambda p: (p / "bin" / "deepcoil").is_file() or (p / "bin" / "python").is_file(),
)

# catGRANULE v1 is reimplemented in src/phasepred/catgranule_v1.py — built-in.

_TOOL_SPECS: list[_ToolSpec] = [
    _SEG, _PSCORE, _PLAAC, _IUPRED3, _ESPRITZ, _DEEPCOIL,
]


# ---------------------------------------------------------------------------
# Public API — one resolver per tool
# ---------------------------------------------------------------------------

def find_seg() -> Path:
    s = _resolve(_SEG)
    if not s.ok:
        raise PhaSePredToolNotFound(f"SEG not found. {s.install_hint}")
    return s.path


def find_pscore_dir() -> Path:
    s = _resolve(_PSCORE)
    if not s.ok:
        raise PhaSePredToolNotFound(f"PScore not found. {s.install_hint}")
    return s.path


def find_pscore_script() -> Path:
    return find_pscore_dir() / "elife_phase_separation_predictor.py"


def find_plaac_wrapper() -> Path:
    """Return the shell wrapper that knows how to run plaac.jar."""
    if not _PLAAC_WRAPPER.is_file():
        raise PhaSePredToolNotFound(
            f"PLAAC wrapper not found at {_PLAAC_WRAPPER}. "
            "Check that the repository is intact under tools/wrappers/."
        )
    # Sanity: also confirm the jar can be resolved so the wrapper won't fail at runtime.
    s = _resolve(_PLAAC)
    if not s.ok:
        raise PhaSePredToolNotFound(f"PLAAC jar not found. {s.install_hint}")
    return _PLAAC_WRAPPER


def find_plaac_jar() -> Path:
    s = _resolve(_PLAAC)
    if not s.ok:
        raise PhaSePredToolNotFound(f"PLAAC jar not found. {s.install_hint}")
    return s.path


def find_iupred3_dir() -> Path:
    s = _resolve(_IUPRED3)
    if not s.ok:
        raise PhaSePredToolNotFound(f"IUPred3 not found. {s.install_hint}")
    return s.path


def find_espritz_wrapper() -> Path:
    """Return the shell wrapper that drives ESpritz from its install dir."""
    wrapper = WRAPPERS_ROOT / "run_espritz.sh"
    if not wrapper.is_file():
        raise PhaSePredToolNotFound(
            f"ESpritz wrapper not found at {wrapper}. "
            "Check that the repository is intact under tools/wrappers/."
        )
    s = _resolve(_ESPRITZ)
    if not s.ok:
        raise PhaSePredToolNotFound(f"ESpritz not found. {s.install_hint}")
    return wrapper


def find_espritz_dir() -> Path:
    s = _resolve(_ESPRITZ)
    if not s.ok:
        raise PhaSePredToolNotFound(f"ESpritz not found. {s.install_hint}")
    return s.path


def find_deepcoil_env() -> Path:
    s = _resolve(_DEEPCOIL)
    if not s.ok:
        raise PhaSePredToolNotFound(f"DeepCoil env not found. {s.install_hint}")
    return s.path


# ---------------------------------------------------------------------------
# Diagnostic helpers used by `phasepred check-tools` and by predict() preflight
# ---------------------------------------------------------------------------

def check_all_tools() -> list[ToolStatus]:
    """Return a status row for every external tool, plus the built-in entry."""

    statuses = [_resolve(spec) for spec in _TOOL_SPECS]
    # catGRANULE v1 is built-in; always OK if the module imports
    try:
        from phasepred import catgranule_v1  # noqa: F401
        cat_status = ToolStatus(
            name="catGRANULE",
            ok=True,
            source="built-in",
            path=Path(__file__).parent / "catgranule_v1.py",
            required_for_predict=True,
            install_hint="(paper-formula reimplementation in src/phasepred/catgranule_v1.py)",
        )
    except Exception as e:  # pragma: no cover
        cat_status = ToolStatus(
            name="catGRANULE",
            ok=False,
            source="missing",
            path=None,
            required_for_predict=True,
            install_hint=f"Built-in module failed to import: {e}",
        )
    statuses.append(cat_status)
    return statuses


def assert_predict_requirements() -> list[ToolStatus]:
    """Raise PhaSePredToolNotFound if any predict-required tool is MISSING.

    Returns the full status list so callers can render it.
    """
    statuses = check_all_tools()
    missing_required = [s for s in statuses if not s.ok and s.required_for_predict]
    if missing_required:
        msg_lines = [
            "phasepred predict aborted — required external tool(s) missing:",
            "",
        ]
        for s in missing_required:
            msg_lines.append(f"  - {s.name}: {s.install_hint}")
        msg_lines.append("")
        msg_lines.append("Run `phasepred check-tools` for the full status table.")
        raise PhaSePredToolNotFound("\n".join(msg_lines))
    return statuses


def format_status_table(statuses: list[ToolStatus]) -> str:
    """Render the status list as a fixed-width table for terminal display."""
    rows = []
    rows.append(("TOOL", "STATUS", "REQUIRED", "SOURCE", "PATH"))
    for s in statuses:
        rows.append((
            s.name,
            "OK" if s.ok else "MISSING",
            "yes" if s.required_for_predict else "no",
            s.source if s.ok else "—",
            str(s.path) if s.path else "—",
        ))
    widths = [max(len(r[i]) for r in rows) for i in range(5)]
    out = []
    for idx, r in enumerate(rows):
        out.append("  ".join(r[i].ljust(widths[i]) for i in range(5)).rstrip())
        if idx == 0:
            out.append("  ".join("-" * w for w in widths))
    n_ok = sum(1 for s in statuses if s.ok)
    n_missing = sum(1 for s in statuses if not s.ok)
    n_missing_required = sum(1 for s in statuses if not s.ok and s.required_for_predict)
    out.append("")
    out.append(f"{n_ok} OK  ·  {n_missing} MISSING  ·  {n_missing_required} required-but-missing")
    return "\n".join(out)
