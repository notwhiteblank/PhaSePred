"""Integrity tests for the E5 golden fixture matrix (Tasks 2 and 3).

These validate the committed fixtures themselves -- INDEX.json coverage,
sequence provenance and the DeepCoil threshold logic -- using only the files in
``tests/fixtures`` and the product feature definitions. No private archive and
no external tool is touched, so the module runs in the ``toolfree`` tier.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pytest

from phasepred.features import DEEPCOIL_THRESHOLD, binary_deepcoil

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures"
GOLDEN = FIXTURES / "golden"
SEQUENCES = FIXTURES / "sequences"

ARCHIVE_SEQ_SOURCE = "archive-2022-02-11"
PSP_MARKER = "042026"

# Task 3: the four-mode frozen prediction matrix and spec 2.2's published model
# sha256 values (Ruling 3's environment fingerprint is compared against these).
PREDICT_MODES = ("SaPS", "PdPS", "hSaPS", "hPdPS")
PUBLISHED_MODEL_SHA256 = {
    "SaPS": "a369e516fbaf92fa3392ada527f3519963038d928c275a4b263755fc20002835",
    "PdPS": "f82eff5d5da9277ec225a596f1594a826d1b9311356e7181827ff887395c8607",
    "hSaPS": "71314225092dcde39e592106272de6afdc6a51371c995e98ab7d0f2d7729a7e0",
    "hPdPS": "1e42b7b714ea3faf0f862b477655de6af985fd048076cdc8da9527bf393166d4",
}

REQUIRED_FEATURES = {
    "Hydropathy",
    "FCR",
    "IDR",
    "LCR",
    "PScore",
    "PLAAC",
    "DeepCoil",
    "catGRANULE",
    "Phos freq",
    "DeepPhase",
}

# spec 4.5's mandatory boundary list (the DeepCoil proximity case is extra).
REQUIRED_BOUNDARIES = {
    "LCR>0",
    "DeepCoil=1",
    "PLAAC_prion_like_positive",
    "IDR_high",
    "IDR=0",
    "PSP_has_phosphosites",
    "PSP_no_record",
    "DeepPhase_has_record",
    "DeepPhase_no_record",
    "very_short_sequence",
    "nonstandard_residue",
}

pytestmark = pytest.mark.toolfree


@pytest.fixture(scope="module")
def index() -> dict:
    return json.loads((GOLDEN / "INDEX.json").read_text(encoding="utf-8"))


def _golden(acc: str) -> dict:
    return json.loads((GOLDEN / f"{acc}.json").read_text(encoding="utf-8"))


def test_feature_matrix_has_at_least_three_cases_each(index: dict) -> None:
    # A bare total is not enough (E1-E4 lesson): each feature must carry >=3
    # distinct accessions, each with a non-empty reason.
    for feature in REQUIRED_FEATURES:
        cases = index["features"][feature]["accessions"]
        assert len(cases) >= 3, f"{feature} has only {len(cases)} cases"
        accessions = [case["accession"] for case in cases]
        assert len(set(accessions)) == len(accessions), f"{feature} repeats a case"
        for case in cases:
            assert case["accession"].strip()
            assert case["reason"].strip(), f"{feature}/{case['accession']} has no reason"


def test_feature_accessions_resolve_to_committed_fixtures(index: dict) -> None:
    for feature in REQUIRED_FEATURES:
        for case in index["features"][feature]["accessions"]:
            acc = case["accession"]
            assert (GOLDEN / f"{acc}.json").is_file(), acc
            assert (SEQUENCES / f"{acc}.fasta").is_file(), acc
            meta = index["proteins"].get(acc)
            assert meta is not None, f"{acc} missing from INDEX proteins"
            assert meta["role"].strip(), f"{acc} has no role"


def test_boundary_coverage_is_complete_and_identifiable(index: dict) -> None:
    coverage = index["boundary_coverage"]
    assert REQUIRED_BOUNDARIES <= set(coverage)
    assert "DeepCoil_near_threshold_0.82" in coverage
    known = set(index["proteins"])
    for name, accessions in coverage.items():
        assert accessions, f"boundary {name} has no accession"
        for acc in accessions:
            assert acc in known, f"boundary {name} points at unknown {acc}"


def test_boundary_coverage_membership_semantics(index: dict) -> None:
    # T2 review: identifying the boundary accessions is not enough -- assert the
    # values that make each boundary true, the way the DeepCoil test already
    # does. A boundary that points at accessions that do not satisfy it would
    # otherwise pass.
    coverage = index["boundary_coverage"]

    def value(acc: str, tool: str, field: str) -> object:
        return _golden(acc)["tools"][tool][field]

    for acc in coverage["IDR=0"]:
        assert value(acc, "espritz", "single") == 0.0, acc
    for acc in coverage["IDR_high"]:
        idr = value(acc, "espritz", "single")
        assert idr is not None and idr >= 0.5, acc
    for acc in coverage["LCR>0"]:
        lcr = value(acc, "seg", "single")
        assert lcr is not None and lcr > 0.0, acc
    for acc in coverage["DeepCoil=1"]:
        assert value(acc, "deepcoil", "binary") == 1.0, acc
    for acc in coverage["PLAAC_prion_like_positive"]:
        plaac = value(acc, "plaac", "single")
        assert plaac is not None and plaac > 0.0, acc
    for acc in coverage["PSP_has_phosphosites"]:
        assert value(acc, "phos", "site_count") > 0, acc
    for acc in coverage["PSP_no_record"]:
        assert value(acc, "phos", "site_count") == 0, acc
    for acc in coverage["DeepPhase_has_record"]:
        assert value(acc, "deepphase", "has_record") is True, acc
    for acc in coverage["DeepPhase_no_record"]:
        assert value(acc, "deepphase", "has_record") is False, acc
    for acc in coverage["nonstandard_residue"]:
        assert value(acc, "hydropathy", "single") is None, acc


def test_phos_freq_cases_are_record_only(index: dict) -> None:
    # Pre-step / Ruling 2: the per-accession Phos value is the web archive's
    # 2020-09-08 record and must never be labelled or read as the marker-gated
    # gold. The gold is the Task 3 local frozen CSV.
    feature = index["features"]["Phos freq"]
    assert "record-only" in feature["source"]
    assert "Task 3" in feature["source"]
    for case in feature["accessions"]:
        assert case.get("record_only") is True, case["accession"]


def test_golden_seq_source_and_psp_marker(index: dict) -> None:
    assert index["psp_version_marker"] == PSP_MARKER
    for acc, meta in index["proteins"].items():
        golden = _golden(acc)
        assert golden["psp_version_marker"] == PSP_MARKER, acc
        assert golden["seq_source"] == meta["seq_source"], acc
        if meta.get("synthetic"):
            assert golden["seq_source"] == "synthetic"
        else:
            assert golden["seq_source"] == ARCHIVE_SEQ_SOURCE, acc


def test_idr_fasta_matches_the_archived_sequence() -> None:
    # T2-G5: the IDR-bearing FASTA must be the archive sequence, not the current
    # UniProt one. The migrated web API fixtures carry M(state) archive 2022-02-11.
    for acc in ("Q08211", "P35637"):
        fasta = (SEQUENCES / f"{acc}.fasta").read_text(encoding="utf-8")
        sequence = "".join(
            line for line in fasta.splitlines() if not line.startswith(">")
        )
        web = json.loads((GOLDEN / f"web_{acc}.json").read_text(encoding="utf-8"))
        assert sequence == web["data"]["sequence"], acc
        assert _golden(acc)["seq_source"] == ARCHIVE_SEQ_SOURCE


def test_deepcoil_threshold_and_binary_agree_with_the_product() -> None:
    assert DEEPCOIL_THRESHOLD == 0.82
    for path in sorted(GOLDEN.glob("[A-Z]*.json")):
        if path.name == "INDEX.json":
            continue
        acc = path.stem
        tools = json.loads(path.read_text(encoding="utf-8"))["tools"]
        if "deepcoil" not in tools:
            continue
        raw_max = tools["deepcoil"]["raw_cc_max"]
        if raw_max is None:
            continue
        assert tools["deepcoil"]["binary"] == binary_deepcoil([raw_max]), acc


def test_deepcoil_threshold_boundary_case_sits_close_to_0_82(index: dict) -> None:
    # The archive must provide proteins on both sides of the 0.82 threshold;
    # if it could not, INDEX.json would have to record the closest instead.
    near = index["boundary_coverage"]["DeepCoil_near_threshold_0.82"]
    distances = [
        abs(_golden(acc)["tools"]["deepcoil"]["raw_cc_max"] - 0.82) for acc in near
    ]
    assert min(distances) <= 0.01, "no protein sits near the 0.82 threshold"
    sides = {
        _golden(acc)["tools"]["deepcoil"]["binary"] for acc in near
    }
    assert sides == {0.0, 1.0}, "the near-threshold case does not straddle 0.82"


def test_fixture_tree_stays_under_2mb() -> None:
    # G5: the committed fixture tree must stay small; anything derived from the
    # private archive beyond these fixtures is forbidden.
    total = sum(p.stat().st_size for p in FIXTURES.rglob("*") if p.is_file())
    assert total < 2 * 1024 * 1024, f"tests/fixtures is {total} bytes"


def test_short_sequence_boundary_is_synthetic_and_below_seven() -> None:
    fasta = (SEQUENCES / "SHORT.fasta").read_text(encoding="utf-8")
    sequence = "".join(line for line in fasta.splitlines() if not line.startswith(">"))
    assert 0 < len(sequence) < 7
    assert _golden("SHORT")["seq_source"] == "synthetic"


# ---------------------------------------------------------------------------
# Task 3: four-mode frozen prediction matrix (gates T3-G2 / T3-G3)
# ---------------------------------------------------------------------------


def test_frozen_predictions_cover_four_modes_and_three_proteins(index: dict) -> None:
    predictions = index["predictions"]
    assert predictions["modes"] == list(PREDICT_MODES)
    accessions = predictions["accessions"]
    assert len(accessions) >= 3, predictions
    for mode in PREDICT_MODES:
        for acc in accessions:
            path = GOLDEN / f"predict_{mode}_{acc}.csv"
            assert path.is_file(), f"missing frozen prediction {path.name}"
            lines = path.read_text(encoding="utf-8").splitlines()
            # exactly a header plus one row, and the row is that accession
            assert len(lines) == 2, path.name
            assert lines[0].split(",")[0] == "UniprotEntry", path.name
            cells = lines[1].split(",")
            assert cells[0] == acc, path.name
            float(cells[1])  # the score must parse


def test_environment_fingerprint_is_complete(index: dict) -> None:
    fingerprint = index["environment_fingerprint"]
    assert isinstance(fingerprint, dict), "Ruling 3 fingerprint is missing"
    for key in (
        "generation_date",
        "generating_command",
        "python",
        "xgboost",
        "pandas",
        "numpy",
        "scikit_learn",
        "model_sha256",
        "psp_version_marker",
    ):
        assert key in fingerprint, f"fingerprint missing {key!r}"
    assert fingerprint["psp_version_marker"] == PSP_MARKER
    # Ruling 8: the xgboost version is recorded as a FACT, not a compatibility
    # gate -- a 3.4.1 regeneration produced byte-identical CSVs, so the value
    # must not be pinned here (pinning it would fail a valid regeneration).
    assert str(fingerprint["xgboost"]).strip()
    assert set(fingerprint["model_sha256"]) == set(PREDICT_MODES)


def test_environment_fingerprint_note_does_not_claim_version_specificity(
    index: dict,
) -> None:
    # Ruling 8 reverse-proof guard: the emitted note must not resurrect the
    # retracted "bit-exact only in this environment" claim.
    note = index["environment_fingerprint_note"]
    assert "bit-exact only in this environment" not in note
    assert "Ruling 8" in note


def test_ruling9_o75146_retains_both_deepcoil_values() -> None:
    # Ruling 9(a): the conflicting archive and local frozen values are both
    # retained, with the disagreement and its cause recorded. Neither may be
    # deleted or "aligned".
    deepcoil = _golden("O75146")["tools"]["deepcoil"]
    assert deepcoil["raw_cc_max"] == 0.8199
    assert deepcoil["binary"] == 0.0
    assert deepcoil["local_frozen_raw_cc_max"] == 0.82
    assert deepcoil["local_frozen_binary"] == 1.0
    disagreement = deepcoil["disagreement"]
    assert disagreement["archive_binary"] == 0.0
    assert disagreement["local_frozen_binary"] == 1.0
    assert "0.990948" in disagreement["cause"]


def test_index_deepcoil_archive_is_record_only(index: dict) -> None:
    # Ruling 9(b): the archive comparison carries no hard gate; the gold is the
    # local frozen value. The measured discrete agreement is recorded.
    feature = index["features"]["DeepCoil"]
    assert feature["archive_record_only"] is True
    assert feature["archive_discrete_agreement"] == 0.990948
    assert "local frozen" in feature["tolerance"]
    for case in feature["accessions"]:
        if case["accession"] in {"O75146", "O95153", "Q08211", "P35637"}:
            assert "local_frozen_binary" in case, case["accession"]
        else:
            assert "local_frozen_binary" not in case, case["accession"]


def test_model_sha256_matches_published_spec(index: dict) -> None:
    recorded = index["environment_fingerprint"]["model_sha256"]
    assert recorded == PUBLISHED_MODEL_SHA256
    # ... and the recorded values are the actual committed model files, so the
    # fingerprint cannot silently describe a different artifact.
    for mode, expected in PUBLISHED_MODEL_SHA256.items():
        model = REPO_ROOT / "src" / "phasepred" / "data" / "models" / mode / "8f_model_0.joblib"
        assert hashlib.sha256(model.read_bytes()).hexdigest() == expected, mode


def _split_batched_predict_csv(text: str) -> dict[str, str]:
    """Split a batched ``predict`` CSV into ``{accession: header+row}``.

    This is the comparison mechanism E6/E9 reuse: a live batched run must
    reproduce each frozen per-protein file byte-for-byte.
    """
    lines = text.splitlines(keepends=True)
    header = lines[0]
    return {line.split(",", 1)[0]: header + line for line in lines[1:]}


def test_prediction_comparison_detects_a_tampered_score() -> None:
    # T3-G2 reverse proof (cheap, no predict run): a live row equal to the frozen
    # row compares equal; changing one score makes the comparison fail.
    acc = "Q08211"
    frozen_text = (GOLDEN / f"predict_hSaPS_{acc}.csv").read_text(encoding="utf-8")
    header, row = frozen_text.splitlines()
    live = f"{header}\n{row}\n"
    assert _split_batched_predict_csv(live)[acc] == frozen_text

    cells = row.split(",")
    cells[1] = "0.0000000000000000"
    tampered = f"{header}\n{','.join(cells)}\n"
    assert _split_batched_predict_csv(live)[acc] != tampered


def _frozen_phos_freq(accession: str) -> float:
    text = (GOLDEN / f"predict_hSaPS_{accession}.csv").read_text(encoding="utf-8")
    rows = list(csv.DictReader(io.StringIO(text)))
    return float(rows[0]["Phos freq"])


# ---------------------------------------------------------------------------
# Ruling 2: Phos freq is compared exactly only under a matching PSP marker
# (gate T3-G4); a mismatch must SKIP -- never FAIL, never PASS.
# ---------------------------------------------------------------------------


def test_phos_freq_gate_asserts_exactly_on_marker_match(psp_env) -> None:
    # Matching marker -> the comparison happens and is exact.
    psp_env.matches(2.0, 2.0, frozen_marker=PSP_MARKER, installed_marker=PSP_MARKER)
    with pytest.raises(AssertionError):
        psp_env.matches(
            2.0, 2.0000000001, frozen_marker=PSP_MARKER, installed_marker=PSP_MARKER
        )


def test_phos_freq_gate_skips_on_marker_mismatch(psp_env) -> None:
    # Mismatching marker -> SKIP (the exception pytest uses to record a skip),
    # not an assertion failure.
    with pytest.raises(pytest.skip.Exception):
        psp_env.matches(2.0, 2.0, frozen_marker=PSP_MARKER, installed_marker="999999")


def test_frozen_phos_freq_matches_local_psp(index: dict, psp_env) -> None:
    # End-to-end: the frozen Phos freq gold is the local (042026) value. Reverse
    # proof for the SKIP branch: run pytest with
    # PHASEPRED_PSP_VERSION_MARKER=999999 and this test skips with the reason.
    if psp_env.installed_marker() is None:
        pytest.skip("local PhosphoSitePlus dataset unavailable")
    sequence = "".join(
        line
        for line in (SEQUENCES / "Q08211.fasta").read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    from phasepred.tools import compute_phos_freq

    actual = compute_phos_freq([{"accession": "Q08211", "sequence": sequence}])["Q08211"]
    psp_env.matches(
        _frozen_phos_freq("Q08211"),
        actual,
        frozen_marker=index["psp_version_marker"],
    )
