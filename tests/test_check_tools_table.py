"""The human-readable ``check-tools`` table must carry the remedy for MISSING rows.

The hint already lives in :class:`~phasepred.tool_paths.ToolStatus.install_hint`
and is already emitted by ``check-tools --json`` under ``hint``; before this
module the table renderer dropped it, so a reader who saw ``MISSING`` on screen
had no way to know what to install. These tests are deliberately toolfree: they
build synthetic ``ToolStatus`` rows, so they run in CI's ``-m toolfree`` tier and
guard the renderer against a silent regression.

Reverse proof: each case pairs a MISSING row with an OK row and asserts the OK
row's hint is *absent* — if the renderer fell back to dumping every hint, the
second assertion fails; if it stopped emitting hints, the first does.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from phasepred.tool_paths import ToolStatus, format_status_table

pytestmark = pytest.mark.toolfree

TOOLS_HINT = "Data (DBS) is gitignored: bash tools/PScore/install.sh --check verifies script + DBS."
PACKAGE_HINT = "Bundled with the phasepred wheel (deepphase_scores.tsv)."


def _missing(hint: str) -> ToolStatus:
    return ToolStatus(
        name="PScore",
        ok=False,
        source="missing",
        path=None,
        required_for_predict=True,
        license="CC-BY-4.0",
        redistributable=True,
        install_hint=hint,
    )


def _ok(hint: str) -> ToolStatus:
    return ToolStatus(
        name="SEG",
        ok=True,
        source="vendored",
        path=Path("/site-packages/phasepred/data/tools/SEG/run"),
        required_for_predict=True,
        license="ncbi-public-domain-convention",
        redistributable=True,
        install_hint=hint,
    )


def test_table_prints_install_hint_for_missing_rows_only() -> None:
    table = format_status_table([_missing(TOOLS_HINT), _ok("Build the vendored SEG binary.")])

    assert "Install hints:" in table
    assert "PScore: " + TOOLS_HINT in table
    # Reverse proof: an OK row's hint is not a remedy and must not be printed;
    # if the renderer dumped every hint (or the whole ToolStatus repr) this fails.
    assert "Build the vendored SEG binary." not in table


def test_table_flags_repo_only_installers_for_wheel_and_sdist_users() -> None:
    table = format_status_table([_missing(TOOLS_HINT)])

    # tools/ ships in neither the wheel nor the sdist, so the hint points at a
    # path a PyPI-only user does not have; the renderer must say so.
    assert "ships only in the GitHub repository" in table
    assert "wheel" in table
    assert "sdist" in table


def test_table_omits_the_repo_note_when_no_installer_names_tools() -> None:
    # A package row's remedy (reinstall the wheel) is not a repository path, so
    # the note would be false; it must be conditional.
    table = format_status_table([_missing(PACKAGE_HINT)])

    assert PACKAGE_HINT in table
    assert "ships only in the GitHub repository" not in table


def test_table_uses_dash_path_for_missing_rows() -> None:
    # PScore/ESpritz resolve the vendored `run` wrapper, so `path` exists while
    # the runtime behind it does not. Printing that path next to MISSING reads
    # as if the tool were installed; DeepCoil/PhosphoSitePlus show `—` instead.
    # Construct the missing row *with* a path so this exercises that case.
    missing_with_path = ToolStatus(
        name="PScore",
        ok=False,
        source="missing",
        path=Path("/site-packages/phasepred/data/tools/PScore/run"),
        required_for_predict=True,
        license="CC-BY-4.0",
        redistributable=True,
        install_hint=TOOLS_HINT,
    )

    table = format_status_table([missing_with_path])

    assert "/site-packages" not in table
    assert "PScore  MISSING" in table


def test_table_has_no_hint_block_when_every_row_is_ok() -> None:
    table = format_status_table([_ok("Build the vendored SEG binary.")])

    assert "Install hints:" not in table
    assert "ships only in the GitHub repository" not in table