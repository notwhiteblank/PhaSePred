"""D40 drift guard: the wheel-bundled contract copies must not drift from the
checkout.

This test used to live in ``tests/test_tool_packages.py`` under that module's
blanket ``pytest.mark.tools``. It reads only repository files (manifests,
``run``/``install.sh`` bodies, READMEs) and needs no external tool, so the E6
review (Minor 3) moved it to the ``toolfree`` tier -- the tier E7's CI selects
with ``pytest -m toolfree`` -- otherwise the guard that protects D40's debt would
never run in CI. A checkout-side edit to a contract's ``run`` or ``install.sh``
would again go undetected, which is the exact failure class D40 closed.

The tier is pinned by ``test_drift_guard_stays_in_the_toolfree_tier`` below so a
future ``pytestmark`` edit cannot silently move it back out.

``tools_root()`` falls back to ``src/phasepred/data/tools/`` (Ruling 6), so a
wheel user's check-tools and installers read the bundled copies. Every compared
file is pinned on BOTH sides, so either copy drifting alone fails.
"""

from __future__ import annotations

import hashlib
import importlib
import tomllib
from pathlib import Path

import pytest

from phasepred.data import bundled_path

pytestmark = pytest.mark.toolfree

REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = REPO_ROOT / "tools"

PACKAGE_DIR = {
    "seg": "SEG",
    "plaac": "PLAAC",
    "pscore": "PScore",
    "espritz": "ESpritz",
    "deepcoil": "DeepCoil",
    "localcider": "LocalCIDER",
}

# Compare the WHOLE parsed manifest dict minus an explicit allow-list, rather
# than an allow-list of fields to compare. A manifest field added to the schema
# later is then covered by default, instead of drifting silently as
# runtime/path_executables/[[data]]/[output]/[install] could. Nothing is allowed
# to differ today -- D40 aligned both trees byte for byte -- so the allow-list is
# empty; it exists as the one explicit place to record a future deliberate
# asymmetry.
_MANIFEST_ALLOWED_DIFFERENCES: frozenset[str] = frozenset()

# The script and README bodies carry the interpreter/install logic that actually
# broke in D40; a checkout-side edit to `run`, `install.sh` or a README left the
# bundled copy silently rotting while every manifest field still matched.
#
# LocalCIDER is the one explicit, commented exemption (Ruling 2): its bundled
# `run`/`install.sh`/`README.md` deliberately drop the checkout's repo-relative
# .venv/.pixi interpreter candidates, because the bundled copy's inferred repo
# root is the package's data directory, not the checkout. Its manifest is still
# compared below; only its script and README bodies are exempt.
_SCRIPT_DRIFT_FILES = ("run", "install.sh")
_SCRIPT_DRIFT_EXEMPT = frozenset({"LocalCIDER"})
_README_DRIFT_EXEMPT = frozenset({"LocalCIDER"})


def _manifest_without_allowed(manifest: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in manifest.items()
        if key not in _MANIFEST_ALLOWED_DIFFERENCES
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_bundled_contracts_do_not_drift_from_checkout() -> None:
    # D40 / G9. Reverse proofs: change the bundled PScore sha256 pin -> the
    # manifest comparison fails; edit the checkout's tools/DeepCoil/run -> the
    # hash comparison fails; edit either side's README.md -> the README
    # comparison fails.
    bundled_tools = bundled_path("data", "tools")
    for key, dirname in PACKAGE_DIR.items():
        repo_dir = TOOLS_ROOT / dirname
        bundled_dir = bundled_tools / dirname
        repo_manifest = repo_dir / "manifest.toml"
        bundled_manifest = bundled_dir / "manifest.toml"
        assert repo_manifest.is_file(), f"checkout contract missing: {dirname}"
        assert bundled_manifest.is_file(), f"bundled contract missing: {dirname}"
        repo_raw = tomllib.loads(repo_manifest.read_text(encoding="utf-8"))
        bundled_raw = tomllib.loads(bundled_manifest.read_text(encoding="utf-8"))
        assert repo_raw.get("name") == key, (dirname, repo_raw.get("name"))
        # Full dict, allow-list subtracted: a new schema field is compared
        # automatically instead of being silently ignored.
        assert _manifest_without_allowed(bundled_raw) == _manifest_without_allowed(
            repo_raw
        ), dirname

        if dirname not in _SCRIPT_DRIFT_EXEMPT:
            for filename in _SCRIPT_DRIFT_FILES:
                repo_script = repo_dir / filename
                bundled_script = bundled_dir / filename
                assert repo_script.is_file(), ("checkout script missing", dirname, filename)
                assert bundled_script.is_file(), (
                    "bundled script missing",
                    dirname,
                    filename,
                )
                assert _sha256(bundled_script) == _sha256(repo_script), (
                    dirname,
                    filename,
                )

        # E6 review Minor 4: the original nine-file divergence included README
        # bodies (the root tools/README.md among them), so stopping at manifests
        # plus run/install.sh left them free to drift.
        if dirname not in _README_DRIFT_EXEMPT:
            repo_readme = repo_dir / "README.md"
            bundled_readme = bundled_dir / "README.md"
            assert repo_readme.is_file(), ("checkout README missing", dirname)
            assert bundled_readme.is_file(), ("bundled README missing", dirname)
            assert _sha256(bundled_readme) == _sha256(repo_readme), (dirname, "README.md")

    # The root tools/README.md was one of the nine divergent files D40 closed;
    # it has a bundled counterpart and must stay byte-aligned too.
    repo_root_readme = TOOLS_ROOT / "README.md"
    bundled_root_readme = bundled_tools / "README.md"
    assert repo_root_readme.is_file(), "checkout tools/README.md missing"
    assert bundled_root_readme.is_file(), "bundled tools/README.md missing"
    assert _sha256(bundled_root_readme) == _sha256(repo_root_readme), "tools/README.md"

    # PScore's download pins are a hard sha256 gate in its installer. The
    # full-dict comparison above already pins bundled == checkout; this guards
    # the checkout side from losing them, which would silently disable the gate.
    repo_pscore = tomllib.loads(
        (TOOLS_ROOT / "PScore" / "manifest.toml").read_text(encoding="utf-8")
    )
    for field in ("download_url", "sha256", "sha256_policy"):
        assert repo_pscore.get(field), f"checkout PScore manifest lost {field}"


def test_drift_guard_stays_in_the_toolfree_tier() -> None:
    # E6 review Minor 3: the module-level marker is what puts the guard in E7's
    # `pytest -m toolfree` CI tier. Reverse proof: change it to
    # `pytest.mark.tools`, or drop it, and this fails -- so the guard cannot
    # silently leave CI again. Inspecting the module attribute (rather than the
    # collected session) keeps the check deterministic when only this file is
    # selected.
    module = importlib.import_module("test_tool_contract_drift")
    marks = module.pytestmark
    if not isinstance(marks, (list, tuple)):
        marks = [marks]
    names = {mark.name for mark in marks}
    assert "toolfree" in names, names
    assert "tools" not in names, names
