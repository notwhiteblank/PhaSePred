"""Three-tier resolution tests for optional external data (E4 Task 1).

Covers ``data.user_data_root`` / ``data.optional_data_path`` and the lazy
``tools.deepphase_table_path`` / ``tools.phosphosite_path`` functions, plus a
Ruling 1 guard that ``data_path`` never leaks into the XDG/home tiers.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from phasepred import tools
from phasepred.data import (
    data_path,
    optional_data_path,
    repo_root,
    user_data_root,
)

pytestmark = pytest.mark.toolfree

# The private, gitignored DeepPhase workbook is absent from clean public
# checkouts (E6-E8 CI). Tests that need the repo-tier copy to exist must skip
# rather than fail; this mirrors tests/test_deepphase_tsv.py's needs_xlsx gate.
_REPO_XLSX = (
    repo_root() / "data" / "raw" / "external" / "deepphase" / "extracted" / "tableS3.xlsx"
)
needs_xlsx = pytest.mark.skipif(
    not _REPO_XLSX.is_file(),
    reason="DeepPhase tableS3.xlsx absent (see docs/DATA_SOURCES.md)",
)


def _make(root: Path, *parts: str) -> Path:
    target = root.joinpath(*parts)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("x")
    return target


def test_user_data_root_prefers_phasepred_override(tmp_path, monkeypatch):
    override = tmp_path / "override"
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    assert user_data_root() == override.resolve()


def test_user_data_root_uses_xdg_when_override_unset(tmp_path, monkeypatch):
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_DATA_HOME", str(xdg))
    assert user_data_root() == (xdg / "phasepred").resolve()


def test_user_data_root_falls_back_to_local_share(tmp_path, monkeypatch):
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    assert user_data_root() == (home / ".local" / "share" / "phasepred").resolve()


@needs_xlsx
def test_optional_data_path_all_tiers_prefers_override(tmp_path, monkeypatch):
    override = tmp_path / "override"
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    monkeypatch.setenv("XDG_DATA_HOME", str(xdg))
    tier1 = _make(override, "deepphase", "extracted", "tableS3.xlsx")
    _make(xdg, "phasepred", "deepphase", "extracted", "tableS3.xlsx")
    # Tier 3 (the real checkout) also holds this file; tier 1 must still win.
    repo_extracted = repo_root() / "data" / "raw" / "external" / "deepphase" / "extracted"
    assert (repo_extracted / "tableS3.xlsx").is_file()
    got = optional_data_path("deepphase", "extracted", "tableS3.xlsx")
    assert got == tier1


def test_optional_data_path_prefers_xdg_over_repo(tmp_path, monkeypatch):
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_DATA_HOME", str(xdg))
    tier2 = _make(xdg, "phasepred", "deepphase", "extracted", "tableS3.xlsx")
    got = optional_data_path("deepphase", "extracted", "tableS3.xlsx")
    repo_xlsx = repo_root() / "data" / "raw" / "external" / "deepphase" / "extracted"
    assert got == tier2
    assert got != repo_xlsx / "tableS3.xlsx"


def test_optional_data_path_missing_returns_preferred_target(tmp_path, monkeypatch):
    override = tmp_path / "override"
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    got = optional_data_path("pscore", "DBS", "PScore_flat.txt")
    assert got == (override / "pscore" / "DBS" / "PScore_flat.txt").resolve()
    assert not got.exists()


def test_optional_data_path_accepts_legacy_layout(tmp_path, monkeypatch):
    override = tmp_path / "override"
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    legacy = _make(
        override,
        "data",
        "raw",
        "external",
        "phosphositeplus",
        "Phosphorylation_site_dataset.gz",
    )
    got = optional_data_path("phosphositeplus", "Phosphorylation_site_dataset.gz")
    assert got == legacy


def test_optional_data_path_resolves_relative_xdg(tmp_path, monkeypatch):
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("XDG_DATA_HOME", "relative-xdg")
    target = _make(
        tmp_path, "relative-xdg", "phasepred", "pscore", "DBS", "PScore_flat.txt"
    )
    got = optional_data_path("pscore", "DBS", "PScore_flat.txt")
    assert got.is_absolute()
    assert got == target.resolve()


def test_optional_data_path_creates_no_directories(tmp_path, monkeypatch):
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    before = sorted(tmp_path.rglob("*"))
    got = optional_data_path("pscore", "DBS", "PScore_flat.txt")
    assert not got.exists()
    assert sorted(tmp_path.rglob("*")) == before


def test_optional_data_path_unknown_component_raises_keyerror(tmp_path, monkeypatch):
    # ``_COMPONENT_LEGACY_NAMES`` is an explicit table: an unknown component must
    # fail loudly rather than guess a path by string concatenation.
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    with pytest.raises(KeyError):
        optional_data_path("not-a-component", "file.txt")


def test_phosphosite_path_reads_env_at_call_time(tmp_path, monkeypatch):
    # ``phasepred.tools`` was imported at module load, before this test sets the
    # environment. A module-level constant would have frozen the checkout path;
    # the lazy function must instead see the newly-set override.
    assert "phasepred.tools" in sys.modules
    override = tmp_path / "late-root"
    target = _make(override, "phosphositeplus", "Phosphorylation_site_dataset.gz")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    got = tools.phosphosite_path()
    assert got == target
    assert got.exists()


@pytest.mark.parametrize("with_xdg", [True, False])
def test_data_path_ignores_xdg_and_stays_in_checkout(tmp_path, monkeypatch, with_xdg):
    # Ruling 1 reverse proof: the new XDG/home tier must never leak into
    # data_root()/data_path(), which serve committed artifacts under the repo.
    # Parameterized over ``XDG_DATA_HOME`` set/unset so that a home-only
    # ``~/.local/share`` fallback (which would only fire when XDG is unset)
    # cannot escape; ``HOME`` is faked and also holds a poisoned layout.
    monkeypatch.delenv("PHASEPRED_DATA_ROOT", raising=False)
    xdg = tmp_path / "xdg"
    (xdg / "data" / "processed").mkdir(parents=True)
    (xdg / "data" / "processed" / "fake.tsv").write_text("fake")
    (xdg / "phasepred" / "data" / "processed").mkdir(parents=True)
    (xdg / "phasepred" / "data" / "processed" / "fake.tsv").write_text("fake")
    home = tmp_path / "home"
    home_share = home / ".local" / "share" / "phasepred" / "data" / "processed"
    home_share.mkdir(parents=True)
    (home_share / "fake.tsv").write_text("fake")
    monkeypatch.setenv("HOME", str(home))
    if with_xdg:
        monkeypatch.setenv("XDG_DATA_HOME", str(xdg))
    else:
        monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    got = data_path("data", "processed")
    assert got == repo_root() / "data" / "processed"
    assert got.is_relative_to(repo_root())
    assert not got.is_relative_to(xdg)
    assert not got.is_relative_to(home)


def test_deepphase_table_path_prefers_bundled_tsv(tmp_path, monkeypatch):
    bundled = _make(tmp_path, "data", "deepphase", "deepphase_scores.tsv")
    monkeypatch.setattr(tools, "bundled_path", lambda *parts: bundled)
    assert tools.deepphase_table_path() == bundled


def test_deepphase_table_path_falls_back_without_bundled_tsv(tmp_path, monkeypatch):
    monkeypatch.setattr(
        tools, "bundled_path", lambda *parts: tmp_path / "absent" / "deepphase_scores.tsv"
    )
    override = tmp_path / "override"
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(override))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    xlsx = _make(override, "deepphase", "extracted", "tableS3.xlsx")
    assert tools.deepphase_table_path() == xlsx


def test_lookup_deepphase_warns_and_returns_empty_when_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(tools, "deepphase_table_path", lambda: tmp_path / "absent.xlsx")
    with pytest.warns(tools.PhaSePredMissingDataWarning):
        assert tools.lookup_deepphase(["Q08211"]) == {}
