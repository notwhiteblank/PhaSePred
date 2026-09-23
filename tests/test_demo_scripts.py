"""The demo entry points at the repository root (E8): train.sh and validate.sh.

These tests are the committed form of the reverse proofs. They assert the
unhappy paths, not just the happy one: a corrupted table must stop the preflight
with exit 2, and a deliberately different artifact must make the comparison exit
1 rather than pass silently.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
TRAIN = REPO_ROOT / "train.sh"
VALIDATE = REPO_ROOT / "validate.sh"
MODELS = REPO_ROOT / "src" / "phasepred" / "data" / "models"
TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")

pytestmark = pytest.mark.toolfree


def run_script(
    script: Path, *args: str, workdir: Path | None = None
) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PHASEPRED_PYTHON"] = sys.executable
    return subprocess.run(
        ["bash", str(script), *args],
        cwd=str(workdir or REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
    )


def minimal_env() -> dict[str, str]:
    return {"PATH": "/usr/bin:/bin", "PHASEPRED_PYTHON": "/nonexistent/python"}


@pytest.mark.parametrize("script", [TRAIN, VALIDATE])
def test_help_works_without_a_provisioned_environment(script: Path):
    proc = subprocess.run(
        ["bash", str(script), "--help"], env=minimal_env(), capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
    assert "Exit codes" in proc.stdout


def test_compare_only_self_is_byte_identical(tmp_path: Path):
    work = tmp_path / "work"
    shutil.copytree(MODELS, work / "models")
    proc = run_script(TRAIN, "--compare-only", "--work-dir", str(work), "--models-dir", str(MODELS))
    assert proc.returncode == 0, proc.stderr
    assert "40/40 byte-identical" in proc.stdout


def test_compare_only_flags_a_differing_artifact(tmp_path: Path):
    work = tmp_path / "work"
    shutil.copytree(MODELS, work / "models")
    victim = work / "models" / "SaPS" / "8f_model_0.joblib"
    victim.write_bytes(victim.read_bytes() + b"x")
    proc = run_script(TRAIN, "--compare-only", "--work-dir", str(work), "--models-dir", str(MODELS))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "39/40 byte-identical" in proc.stdout
    assert "SaPS/8f_model_0.joblib" in proc.stdout


def test_smoke_mode_keeps_the_shipped_models_unchanged(tmp_path: Path):
    # Exercises the shipped-artifact snapshot guard on the fast (smoke) path.
    before = {p: p.read_bytes() for p in MODELS.rglob("*.joblib")}
    proc = run_script(TRAIN, "--smoke", "20", "--work-dir", str(tmp_path / "smoke"))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "smoke OK" in proc.stdout
    after = {p: p.read_bytes() for p in MODELS.rglob("*.joblib")}
    assert before == after


def _fake_repo(tmp_path: Path, script_name: str, mutate) -> Path:
    """A minimal repo-like tree whose only defect is what ``mutate`` injects."""
    fake = tmp_path / "fake"
    (fake / "data" / "processed").mkdir(parents=True)
    for table in (REPO_ROOT / "data" / "processed").glob("*.tsv"):
        shutil.copy2(table, fake / "data" / "processed" / table.name)
    (fake / "src").symlink_to(REPO_ROOT / "src", target_is_directory=True)
    shutil.copy2(REPO_ROOT / script_name, fake / script_name)
    mutate(fake)
    return fake


@pytest.mark.parametrize("script_name", ["train.sh", "validate.sh"])
def test_preflight_rejects_a_corrupted_training_table(tmp_path: Path, script_name: str):
    def corrupt(fake: Path) -> None:
        table = fake / "data" / "processed" / "chen2022_s2s3_human-features_test.tsv"
        table.write_bytes(table.read_bytes() + b"x")

    fake = _fake_repo(tmp_path, script_name, corrupt)
    proc = run_script(fake / script_name, "--work-dir", str(tmp_path / "noop"))
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "sha256 mismatch" in proc.stderr


def test_validate_propagates_a_hard_consistency_failure(tmp_path: Path):
    models = tmp_path / "models"
    shutil.copytree(MODELS, models)
    manifest = models / "SaPS" / "manifest.json"
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["feature_definitions_version"] = "v1999"
    manifest.write_text(json.dumps(payload), encoding="utf-8")

    proc = run_script(VALIDATE, "--models-dir", str(models), "--work-dir", str(tmp_path / "work"))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "hard_fail=2" in proc.stdout or "hard_fail=1" in proc.stdout


def test_scripts_have_no_absolute_machine_paths():
    for script in (TRAIN, VALIDATE):
        text = script.read_text(encoding="utf-8")
        assert "/mnt/" not in text
        assert "/home/" not in text
