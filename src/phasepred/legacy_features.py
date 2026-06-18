from __future__ import annotations

import gzip
import json
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

from phasepred.features import compute_native_features
from phasepred.tool_paths import (
    PhaSePredToolNotFound,
    find_plaac_wrapper,
    find_pscore_dir,
    find_pscore_script,
    find_seg,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


class LegacyFeatureError(RuntimeError):
    """Raised when a legacy feature adapter cannot compute a value."""


@dataclass(frozen=True)
class ScalarFeature:
    accession: str
    value: float
    source: str

    def to_json(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_json(cls, payload: dict[str, object]) -> ScalarFeature:
        return cls(
            accession=str(payload["accession"]),
            value=float(payload["value"]),
            source=str(payload["source"]),
        )


class ScalarFeatureCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def read_all(self) -> dict[str, ScalarFeature]:
        if not self.path.exists():
            return {}
        rows: dict[str, ScalarFeature] = {}
        with self.path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if stripped:
                    feature = ScalarFeature.from_json(json.loads(stripped))
                    rows[feature.accession] = feature
        return rows

    def append_many(self, features: Iterable[ScalarFeature]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            for feature in features:
                handle.write(json.dumps(feature.to_json(), sort_keys=True) + "\n")


def compute_lcr_features(sequence_rows: list[dict[str, object]]) -> list[ScalarFeature]:
    try:
        seg_bin = find_seg()
    except PhaSePredToolNotFound as e:
        raise LegacyFeatureError(str(e)) from e
    fasta = _write_temp_fasta(sequence_rows)
    try:
        result = subprocess.run(
            [str(seg_bin), str(fasta), "12", "2.2", "2.5", "-x"],
            check=True,
            capture_output=True,
            text=True,
        )
    finally:
        fasta.unlink(missing_ok=True)
    masked = _parse_fasta_text(result.stdout)
    features = []
    for row in sequence_rows:
        accession = str(row["UniprotEntry"])
        sequence = masked.get(accession, "")
        if not sequence:
            continue
        low_complexity = sum(1 for residue in sequence if residue == "x" or residue.islower())
        features.append(
            ScalarFeature(
                accession=accession,
                value=low_complexity / len(sequence),
                source="seg_12_2.2_2.5_x",
            )
        )
    return features


def compute_pscore_features(sequence_rows: list[dict[str, object]]) -> list[ScalarFeature]:
    try:
        pscore_script = find_pscore_script()
        pscore_dir = find_pscore_dir()
    except PhaSePredToolNotFound as e:
        raise LegacyFeatureError(str(e)) from e
    with tempfile.TemporaryDirectory(prefix="phasepred_pscore_") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "input.fasta"
        output = workdir / "pscore.tsv"
        fasta.write_text(_fasta_text(sequence_rows), encoding="utf-8")
        subprocess.run(
            [
                "uv",
                "run",
                "python",
                str(pscore_script),
                str(fasta),
                "-output",
                str(output),
                "-overwrite",
                "-mute",
            ],
            cwd=pscore_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        return parse_pscore_output(output)


def parse_pscore_output(path: str | Path) -> list[ScalarFeature]:
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].startswith("PScore"):
            accession = parts[-1].lstrip(">")
            rows.append(ScalarFeature(accession, float(parts[1]), "pscore_elife_sourcecode2"))
    return rows


def compute_plaac_features(sequence_rows: list[dict[str, object]]) -> list[ScalarFeature]:
    try:
        plaac_runner = find_plaac_wrapper()
    except PhaSePredToolNotFound as e:
        raise LegacyFeatureError(str(e)) from e
    with tempfile.TemporaryDirectory(prefix="phasepred_plaac_") as tmp:
        fasta = Path(tmp) / "input.fasta"
        fasta.write_text(_fasta_text(sequence_rows), encoding="utf-8")
        result = subprocess.run(
            [str(plaac_runner), "-i", str(fasta)],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    return parse_plaac_output(result.stdout)


def parse_plaac_output(text: str) -> list[ScalarFeature]:
    lines = [line for line in text.splitlines() if line.strip()]
    header = None
    rows = []
    for line in lines:
        if line.startswith("SEQid\t"):
            header = line.split("\t")
            continue
        if header is None or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != len(header):
            continue
        record = dict(zip(header, parts, strict=False))
        if "SEQid" in record and "NLLR" in record:
            rows.append(ScalarFeature(record["SEQid"], float(record["NLLR"]), "plaac_nllr"))
    return rows


def compute_phosphosite_frequencies(
    sequence_rows: list[dict[str, object]],
    phosphosite_path: str | Path,
) -> list[ScalarFeature]:
    lengths = {str(row["UniprotEntry"]): len(str(row["sequence"])) for row in sequence_rows}
    counts = dict.fromkeys(lengths, 0)
    with gzip.open(phosphosite_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("GENE\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None or not line.strip() or line.startswith("PhosphoSitePlus"):
                continue
            record = dict(zip(header, line.rstrip("\n").split("\t"), strict=False))
            accession = record.get("ACC_ID", "")
            if accession in counts and record.get("ORGANISM", "").lower() == "human":
                counts[accession] += 1
    return [
        ScalarFeature(accession, count / lengths[accession], "phosphositeplus_sites_per_residue")
        for accession, count in counts.items()
        if lengths[accession] > 0
    ]


def load_deepphase_scores(path: str | Path) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name="tableS3")
    required = {"Swiss-Prot ID", "DeepPhase_score"}
    missing = required - set(frame.columns)
    if missing:
        raise LegacyFeatureError(f"DeepPhase table missing columns: {', '.join(sorted(missing))}")
    return frame[["Swiss-Prot ID", "DeepPhase_score"]].rename(
        columns={"Swiss-Prot ID": "UniprotEntry", "DeepPhase_score": "DeepPhase"}
    )


def native_feature_frame(sequence_rows: list[dict[str, object]]) -> pd.DataFrame:
    rows = []
    errors = []
    for row in sequence_rows:
        accession = str(row["UniprotEntry"])
        try:
            features = compute_native_features(str(row["sequence"]))
        except ValueError as exc:
            errors.append({"UniprotEntry": accession, "error": str(exc)})
            features = {}
        rows.append({"UniprotEntry": accession, **features})
    frame = pd.DataFrame(rows)
    if errors:
        frame.attrs["errors"] = errors
    return frame


def _write_temp_fasta(sequence_rows: list[dict[str, object]]) -> Path:
    handle = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8")
    with handle:
        handle.write(_fasta_text(sequence_rows))
    return Path(handle.name)


def _fasta_text(sequence_rows: list[dict[str, object]]) -> str:
    return "".join(
        f">{str(row['UniprotEntry']).strip()}\n{''.join(str(row['sequence']).split()).upper()}\n"
        for row in sequence_rows
    )


def _parse_fasta_text(text: str) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return {accession: "".join(parts) for accession, parts in records.items()}
