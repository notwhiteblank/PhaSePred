"""Driver-level tests for scripts/train_phasepred.py (E2).

Loads the driver the same way tests/test_build_training_tsv.py does and covers
the smoke-target guard: smoke mode truncates the negative pool, so its metrics
and accession lists are not the published ones and must never land on the
published artifact locations. None of these tests train anything.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "train_phasepred.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("train_phasepred", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["train_phasepred"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def driver():
    return load_driver()


def test_resolve_defaults_unchanged(driver):
    models, metrics = driver.resolve_output_dirs(None, None, None)
    assert models == driver.models_root()
    assert metrics == driver.data_path("products", "A_paper_split_recomputed")


def test_smoke_refuses_without_a_metrics_dir(driver, tmp_path):
    with pytest.raises(SystemExit) as exc:
        driver.resolve_output_dirs(tmp_path, None, 10)
    assert exc.value.code == 2


def test_smoke_refuses_without_an_output_dir(driver, tmp_path):
    with pytest.raises(SystemExit) as exc:
        driver.resolve_output_dirs(None, tmp_path, 10)
    assert exc.value.code == 2


def test_smoke_accepts_two_temporary_dirs(driver, tmp_path):
    models = tmp_path / "models"
    metrics = tmp_path / "metrics"
    assert driver.resolve_output_dirs(models, metrics, 10) == (models, metrics)


def test_smoke_refuses_the_published_product_dir(driver, tmp_path):
    published = driver.data_path("products", "A_paper_split_recomputed")
    with pytest.raises(SystemExit):
        driver.resolve_output_dirs(tmp_path, published, 10)


def test_smoke_refuses_the_published_models_dir(driver, tmp_path):
    with pytest.raises(SystemExit):
        driver.resolve_output_dirs(driver.models_root(), tmp_path, 10)


def test_smoke_refuses_a_subdirectory_of_the_published_models_dir(driver, tmp_path):
    with pytest.raises(SystemExit):
        driver.resolve_output_dirs(driver.models_root() / "SaPS", tmp_path, 10)


def test_smoke_refuses_metrics_inside_the_published_models_dir(driver, tmp_path):
    # Cross case: a temporary --output-dir with --metrics-dir landing inside the
    # published models root would overwrite a per-task manifest from smoke data.
    with pytest.raises(SystemExit) as exc:
        driver.resolve_output_dirs(tmp_path, driver.models_root() / "SaPS", 10)
    assert exc.value.code == 2


def test_smoke_refuses_output_inside_the_published_models_dir(driver, tmp_path):
    with pytest.raises(SystemExit) as exc:
        driver.resolve_output_dirs(driver.models_root() / "SaPS", tmp_path, 10)
    assert exc.value.code == 2


def test_smoke_accepts_two_distinct_temporary_dirs(driver, tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    assert driver.resolve_output_dirs(first, second, 10) == (first, second)


def test_refusal_message_names_the_offending_path(driver, tmp_path, capsys):
    published = driver.data_path("products", "A_paper_split_recomputed")
    with pytest.raises(SystemExit):
        driver.resolve_output_dirs(tmp_path, published, 10)
    err = capsys.readouterr().err
    assert str(published) in err
    assert "truncated" in err
