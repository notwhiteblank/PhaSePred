from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol

import pandas as pd
import requests
from Bio import SeqIO

UNIPROT_FASTA_BASE_URL = "https://rest.uniprot.org/uniprotkb"


@dataclass(frozen=True)
class ProteinSequence:
    accession: str
    sequence: str
    source: str
    header: str
    retrieved_at: str | None = None
    endpoint: str | None = None
    sequence_hash: str = field(init=False)
    length: int = field(init=False)

    def __post_init__(self) -> None:
        normalized = self.sequence.replace(" ", "").replace("\n", "").upper()
        if not normalized:
            raise ValueError(f"Protein sequence for {self.accession!r} is empty")
        object.__setattr__(self, "sequence", normalized)
        object.__setattr__(self, "length", len(normalized))
        object.__setattr__(
            self,
            "sequence_hash",
            hashlib.sha256(normalized.encode("utf-8")).hexdigest(),
        )

    def to_json(self) -> dict[str, object]:
        payload = asdict(self)
        payload.pop("sequence_hash", None)
        payload.pop("length", None)
        payload["sequence_hash"] = self.sequence_hash
        payload["length"] = self.length
        return payload

    @classmethod
    def from_json(cls, payload: dict[str, object]) -> ProteinSequence:
        return cls(
            accession=str(payload["accession"]),
            sequence=str(payload["sequence"]),
            source=str(payload["source"]),
            header=str(payload.get("header") or ""),
            retrieved_at=payload.get("retrieved_at") and str(payload["retrieved_at"]),
            endpoint=payload.get("endpoint") and str(payload["endpoint"]),
        )


class SequenceFetcher(Protocol):
    def fetch_fasta(self, accession: str) -> ProteinSequence: ...


class SequenceCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def get(self, accession: str) -> ProteinSequence | None:
        accession = accession.strip()
        for record in self._read_all().values():
            if record.accession == accession:
                return record
        return None

    def put(self, record: ProteinSequence) -> None:
        records = self._read_all()
        records[record.accession] = record
        self.write_all(records)

    def write_all(self, records: dict[str, ProteinSequence]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("w", encoding="utf-8") as handle:
            for item in records.values():
                handle.write(json.dumps(item.to_json(), sort_keys=True) + "\n")

    def _read_all(self) -> dict[str, ProteinSequence]:
        if not self.path.exists():
            return {}
        records: dict[str, ProteinSequence] = {}
        with self.path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                record = ProteinSequence.from_json(json.loads(stripped))
                records[record.accession] = record
        return records


class UniProtClient:
    def __init__(
        self,
        *,
        base_url: str = UNIPROT_FASTA_BASE_URL,
        session: requests.Session | None = None,
        timeout: float = 30.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()
        self.timeout = timeout

    def fetch_fasta(self, accession: str) -> ProteinSequence:
        accession = accession.strip()
        endpoint = f"{self.base_url}/{accession}.fasta"
        response = self.session.get(
            endpoint,
            headers={"Accept": "text/x-fasta"},
            timeout=self.timeout,
        )
        response.raise_for_status()
        records = _records_from_fasta_text(response.text, source="uniprot", endpoint=endpoint)
        if not records:
            raise ValueError(f"UniProt returned no FASTA record for {accession}")
        record = records[0]
        return ProteinSequence(
            accession=accession,
            sequence=record.sequence,
            source="uniprot",
            header=record.header,
            retrieved_at=datetime.now(UTC).isoformat(),
            endpoint=endpoint,
        )


def read_fasta_records(path: str | Path) -> list[ProteinSequence]:
    text = Path(path).read_text(encoding="utf-8")
    return _records_from_fasta_text(text, source="fasta", endpoint=None)


def resolve_accessions(
    accessions: list[str],
    *,
    cache_path: str | Path,
    client: SequenceFetcher | None = None,
    allow_network: bool = True,
) -> list[ProteinSequence]:
    cache = SequenceCache(cache_path)
    fetcher = client or UniProtClient()
    records: list[ProteinSequence] = []
    for accession in accessions:
        normalized = accession.strip()
        if not normalized:
            continue
        cached = cache.get(normalized)
        if cached is not None:
            records.append(cached)
            continue
        if not allow_network:
            raise ValueError(f"Sequence for {normalized} is not cached and network is disabled")
        fetched = fetcher.fetch_fasta(normalized)
        cache.put(fetched)
        records.append(fetched)
    return records


def sequence_records_to_frame(records: list[ProteinSequence]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "UniprotEntry": record.accession,
                "sequence": record.sequence,
                "source": record.source,
                "header": record.header,
                "retrieved_at": record.retrieved_at,
                "endpoint": record.endpoint,
                "length": record.length,
                "sequence_hash": record.sequence_hash,
            }
            for record in records
        ]
    )


def _records_from_fasta_text(
    text: str,
    *,
    source: str,
    endpoint: str | None,
) -> list[ProteinSequence]:
    from io import StringIO

    records: list[ProteinSequence] = []
    for seq_record in SeqIO.parse(StringIO(text), "fasta"):
        header = f">{seq_record.description}"
        accession = _accession_from_header(seq_record.id)
        records.append(
            ProteinSequence(
                accession=accession,
                sequence=str(seq_record.seq),
                source=source,
                header=header,
                endpoint=endpoint,
            )
        )
    if not records:
        raise ValueError("No FASTA records found")
    return records


def _accession_from_header(identifier: str) -> str:
    parts = identifier.split("|")
    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
        return parts[1]
    match = re.match(r"([A-Za-z0-9_.-]+)", identifier)
    if not match:
        raise ValueError(f"Cannot parse FASTA accession from header {identifier!r}")
    return match.group(1)
