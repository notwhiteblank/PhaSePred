from __future__ import annotations

from pathlib import Path

import pytest

from phasepred.sequences import (
    ProteinSequence,
    SequenceCache,
    UniProtClient,
    read_fasta_records,
    resolve_accessions,
)


def test_read_fasta_records_parse_accession_and_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">sp|P12345|GENE_HUMAN example\nACDEK\n", encoding="utf-8")

    records = read_fasta_records(fasta)

    assert len(records) == 1
    record = records[0]
    assert record.accession == "P12345"
    assert record.sequence == "ACDEK"
    assert record.length == 5
    assert record.source == "fasta"


def test_sequence_cache_round_trip(tmp_path: Path) -> None:
    cache = SequenceCache(tmp_path / "cache.jsonl")
    record = ProteinSequence(
        accession="P12345",
        sequence="ACDEK",
        source="fasta",
        header=">sp|P12345|GENE_HUMAN example",
    )

    cache.put(record)
    loaded = SequenceCache(tmp_path / "cache.jsonl").get("P12345")

    assert loaded is not None
    assert loaded.accession == "P12345"
    assert loaded.sequence == "ACDEK"


def test_uniprot_client_fetches_and_parses_fasta_without_network() -> None:
    class DummyResponse:
        status_code = 200
        text = ">sp|P12345|GENE_HUMAN example\nACDEK\n"

        def raise_for_status(self) -> None:
            return None

    class DummySession:
        def __init__(self) -> None:
            self.calls: list[tuple[str, dict[str, str]]] = []

        def get(self, url: str, headers: dict[str, str], timeout: float) -> DummyResponse:
            self.calls.append((url, headers))
            assert timeout == 30.0
            return DummyResponse()

    session = DummySession()
    client = UniProtClient(session=session)
    record = client.fetch_fasta("P12345")

    assert record.accession == "P12345"
    assert record.source == "uniprot"
    assert session.calls[0][0].endswith("/uniprotkb/P12345.fasta")
    assert session.calls[0][1]["Accept"] == "text/x-fasta"


def test_resolve_accessions_uses_cache_before_client(tmp_path: Path) -> None:
    cache_path = tmp_path / "cache.jsonl"
    cache = SequenceCache(cache_path)
    cache.put(
        ProteinSequence(
            accession="CACHED",
            sequence="ACDE",
            source="fasta",
            header=">cached",
        )
    )

    class ExplodingClient:
        def fetch_fasta(self, accession: str) -> ProteinSequence:
            raise AssertionError(f"unexpected network fetch for {accession}")

    records = resolve_accessions(["CACHED"], cache_path=cache_path, client=ExplodingClient())

    assert [record.accession for record in records] == ["CACHED"]


def test_read_fasta_records_rejects_empty_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "bad.fasta"
    fasta.write_text(">sp|P12345|GENE_HUMAN example\n", encoding="utf-8")

    with pytest.raises(ValueError):
        read_fasta_records(fasta)
