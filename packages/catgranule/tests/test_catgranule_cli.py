from __future__ import annotations

import json

from catgranule.cli import app
from typer.testing import CliRunner

runner = CliRunner()


def test_cli_help_hits_package_entrypoint() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "Usage: " in result.output
    assert "predict" in result.output.lower()


def test_cli_predict_json(tmp_path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\nMAPLLLL\n", encoding="utf-8")

    result = runner.invoke(app, ["predict", "--fasta", str(fasta), "--out", "json"])

    assert result.exit_code == 0, result.output
    rows = json.loads(result.output)
    assert len(rows) == 1
    assert rows[0]["accession"] == "P1"
    assert set(rows[0]) == {"accession", "single", "residue", "granule_strength"}
    assert isinstance(rows[0]["single"], float)
    assert isinstance(rows[0]["residue"], list)


def test_cli_predict_json_reports_input_errors() -> None:
    import tempfile
    from pathlib import Path

    from typer.testing import CliRunner as CR

    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "in.fasta"
        fasta.write_text(">P1\nACX\n", encoding="utf-8")
        result = CR().invoke(app, ["predict", "--fasta", str(fasta), "--out", "json"])

    assert result.exit_code == 0
    rows = json.loads(result.output)
    assert rows[0]["single"] is None
    assert "error" in rows[0]


def test_cli_predict_csv(tmp_path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\nMAPLLLL\n", encoding="utf-8")

    result = runner.invoke(app, ["predict", "--fasta", str(fasta), "--out", "csv"])

    assert result.exit_code == 0, result.output
    assert result.output.splitlines()[0] == "accession,single,granule_strength,residue_length"
    assert result.output.splitlines()[1].startswith("P1,")


def test_cli_predict_raw_flag(tmp_path) -> None:
    # long enough sequence so the 51-window interior profile is non-empty
    sequence = "MGGGGGSSSSQQQQKKKKRRRRPPPPHHHHYYYY" + "MGGGGGSSSS" * 7  # 92 res
    fasta = tmp_path / "in.fasta"
    fasta.write_text(f">P1\n{sequence}\n", encoding="utf-8")

    normal_invoke = runner.invoke(app, ["predict", "--fasta", str(fasta), "--out", "json"])
    raw_invoke = runner.invoke(
        app, ["predict", "--fasta", str(fasta), "--out", "json", "--raw"]
    )
    normal = json.loads(normal_invoke.output)
    raw = json.loads(raw_invoke.output)

    # distilled (default) route: --raw keeps the surrogate single, un-aligns
    # the residue profile
    assert normal[0]["single"] == raw[0]["single"]
    assert normal[0]["residue"] != raw[0]["residue"]


def test_cli_predict_route_paper_raw_flag_differs_from_normalized(tmp_path) -> None:
    sequence = "MGGGGGSSSSQQQQKKKKRRRRPPPPHHHHYYYY" + "MGGGGGSSSS" * 7
    fasta = tmp_path / "in.fasta"
    fasta.write_text(f">P1\n{sequence}\n", encoding="utf-8")

    normal = runner.invoke(app, ["predict", "--fasta", str(fasta), "--route", "paper"]).output
    raw = runner.invoke(
        app, ["predict", "--fasta", str(fasta), "--route", "paper", "--raw"]
    ).output
    # paper (formula) route keeps the raw-vs-normalized distinction on single
    assert json.loads(normal)[0]["single"] != json.loads(raw)[0]["single"]


def test_cli_predict_route_switches_artifacts(tmp_path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\nMAPLLLL\n", encoding="utf-8")

    dist = runner.invoke(app, ["predict", "--fasta", str(fasta), "--out", "json"])
    paper = runner.invoke(
        app, ["predict", "--fasta", str(fasta), "--out", "json", "--route", "paper"]
    )
    legacy = runner.invoke(
        app, ["predict", "--fasta", str(fasta), "--out", "json", "--route", "legacy-unaudited"]
    )
    assert dist.exit_code == 0 and paper.exit_code == 0 and legacy.exit_code == 0
    d, p, leg = json.loads(dist.output), json.loads(paper.output), json.loads(legacy.output)
    assert d[0]["single"] != p[0]["single"]
    assert p[0]["single"] != leg[0]["single"]


def test_cli_predict_rejects_unknown_route(tmp_path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">P1\nMAPLLLL\n", encoding="utf-8")
    result = runner.invoke(app, ["predict", "--fasta", str(fasta), "--route", "bogus"])
    assert result.exit_code != 0


def test_cli_check_reports_default_distilled(tmp_path) -> None:
    result = runner.invoke(app, ["check"])
    assert result.exit_code == 0, result.output
    assert "distilled" in result.output
    assert "mode: distilled" in result.output


def test_cli_check_distinguishes_artifact_and_package_versions() -> None:
    # Finding 3 (D43): the distilled artifact stamps 1.0.0 (it was produced by
    # catgranule 1.0.0 and has not changed since; 1.0.1 was a declaration-only
    # convergence), so the label must not read as a package/artifact mismatch.
    import catgranule

    result = runner.invoke(app, ["check"])
    assert result.exit_code == 0, result.output
    assert "artifact_version: 1.0.0" in result.output
    assert f"package_version: {catgranule.__version__}" in result.output
    assert catgranule.__version__ != "1.0.0"