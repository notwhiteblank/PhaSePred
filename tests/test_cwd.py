"""CWD-independence tests (S2): the package must run from any directory.

Two levels are verified:
1. API level (fast, no external tools): data-root resolution and the
   XGBoost predict path must not depend on ``os.getcwd()``.
2. CLI level (subprocess from a foreign CWD): ``phasepred
   features-from-fasta`` through the installed console script, which
   exercises the tool runners (ESpritz + SEG) resolved via absolute paths.

The CLI subprocess runs through :func:`conftest.run_cli`; the runner is the
installed console script (Ruling 5), overridable with ``PHASEPRED_RUNNER``.
``python -m phasepred.cli`` is a silent no-op and is deliberately not used —
``test_cli_runner_produces_parseable_output`` fails if the runner ever
regresses to a form that emits nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from conftest import run_cli

from phasepred.data import data_path, data_root, models_root

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.toolfree
def test_data_root_resolution_is_cwd_independent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)

    assert data_root() == REPO_ROOT
    assert (models_root() / "SaPS").is_dir()
    assert data_path("pyproject.toml").exists()


@pytest.mark.toolfree
def test_data_root_honors_phasespred_data_root_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PHASEPRED_DATA_ROOT", str(tmp_path))

    assert data_root() == tmp_path
    assert data_path("somefile.txt") == tmp_path / "somefile.txt"


@pytest.mark.toolfree
def test_predict_scores_from_foreign_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Full predict path (model load + XGB scoring) must run from /tmp."""
    import numpy as np

    import phasepred.predictor as predictor
    from phasepred.features import BASE_FEATURE_COLUMNS

    monkeypatch.chdir(tmp_path)

    model_dir = models_root() / "SaPS"
    models = predictor.load_ensemble(model_dir)
    features = pd.DataFrame(
        [
            {
                "UniprotEntry": "Q08211",
                "length": 1270.0,
                "Hydropathy": 0.4623,
                "FCR": 0.2252,
                "IDR": 0.0,
                "LCR": 0.1126,
                "PScore": 8.96,
                "PLAAC": 0.224,
                "catGRANULE": 1.5207,
                "DeepCoil": 0.0,
            }
        ]
    )

    predictions = predictor.predict_scores(models, features, list(BASE_FEATURE_COLUMNS))

    assert list(predictions.columns) == ["UniprotEntry", "score"]
    assert np.isfinite(predictions["score"][0])


@pytest.mark.tools
def test_cli_features_from_fasta_from_foreign_cwd(tmp_path: Path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\n" + "M" + "S" * 150 + "\n", encoding="utf-8")
    output = tmp_path / "out.csv"

    proc = run_cli(
        "features-from-fasta",
        "--input",
        str(fasta),
        "--output",
        str(output),
        cwd=tmp_path,
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip(), "CLI produced no stdout (silent-no-op runner?)"
    frame = pd.read_csv(output)
    assert list(frame["UniprotEntry"]) == ["P1"]
    assert "espritz-idr-fraction" in frame.columns
    assert "lcr-fraction" in frame.columns


@pytest.mark.toolfree
def test_cli_runner_produces_parseable_output() -> None:
    # Ruling 5 / T1-G3: the runner must be a form that actually runs. `python
    # -m phasepred.cli` exits 0 with EMPTY stdout, so a returncode-only
    # assertion would pass vacuously. Requiring non-empty, JSON-parseable
    # output fails loudly if the default ever regresses to that no-op.
    # Reverse proof: PHASEPRED_RUNNER="<python> -m phasepred.cli" -> empty
    # stdout -> json.loads raises -> this test fails.
    proc = run_cli("check-tools", "--json")

    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert isinstance(payload, list)
    assert len(payload) == 9
