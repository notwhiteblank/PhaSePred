"""Version declarations must agree (E2, Ruling 1)."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest
from conftest import run_cli

import phasepred

REPO_ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.toolfree


def read_version(pyproject: Path) -> str:
    with pyproject.open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def test_root_version_agrees_everywhere():
    declared = read_version(REPO_ROOT / "pyproject.toml")
    assert phasepred.__version__ == declared
    pixi = (REPO_ROOT / "pixi.toml").read_text(encoding="utf-8")
    assert re.search(rf'^version = "{re.escape(declared)}"$', pixi, re.M), pixi[:200]


def test_cli_version_prints_the_declared_version():
    # Finding 5 (D43): `phasepred --version` must report the same string as
    # ``phasepred.__version__`` and ``pyproject.toml``, deriving from the former
    # so it is not a fourth independent version source.
    declared = read_version(REPO_ROOT / "pyproject.toml")
    proc = run_cli("--version")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == declared
    assert proc.stdout.strip() == phasepred.__version__


def test_catgranule_version_agrees_everywhere():
    import catgranule

    declared = read_version(REPO_ROOT / "packages" / "catgranule" / "pyproject.toml")
    assert catgranule.__version__ == declared


def test_requires_python_admits_313():
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        requires = tomllib.load(handle)["project"]["requires-python"]
    assert requires == ">=3.12"


def test_both_packages_admit_the_same_python_versions():
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        root = tomllib.load(handle)["project"]["requires-python"]
    with (REPO_ROOT / "packages" / "catgranule" / "pyproject.toml").open("rb") as handle:
        sub = tomllib.load(handle)["project"]["requires-python"]
    assert root == sub, f"phasepred {root!r} vs catgranule {sub!r}"


def test_xgboost_pin_excludes_the_non_reproducing_minor():
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        deps = tomllib.load(handle)["project"]["dependencies"]
    xgb = [d for d in deps if d.startswith("xgboost")]
    assert xgb == ["xgboost>=3.2,<3.3"], xgb
