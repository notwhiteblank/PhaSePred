from __future__ import annotations

import json
import subprocess
import tempfile
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

ESPRITZ_RUNNER = Path(__file__).resolve().parents[2] / "scripts" / "tools" / "run_espritz.sh"


@dataclass(frozen=True)
class ESpritzFeature:
    accession: str
    espritz_idr: float
    mean_score: float
    max_score: float
    disorder_threshold: float
    residues: int
    model: str
    sw: int

    def to_json(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_json(cls, payload: dict[str, object]) -> ESpritzFeature:
        return cls(
            accession=str(payload["accession"]),
            espritz_idr=float(payload["espritz_idr"]),
            mean_score=float(payload["mean_score"]),
            max_score=float(payload["max_score"]),
            disorder_threshold=float(payload["disorder_threshold"]),
            residues=int(payload["residues"]),
            model=str(payload["model"]),
            sw=int(payload["sw"]),
        )


class ESpritzFeatureCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def _read_all(self) -> dict[str, ESpritzFeature]:
        if not self.path.exists():
            return {}
        features: dict[str, ESpritzFeature] = {}
        with self.path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                feature = ESpritzFeature.from_json(json.loads(stripped))
                features[feature.accession] = feature
        return features

    def append_new(self, feature: ESpritzFeature) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(feature.to_json(), sort_keys=True) + "\n")


def resolve_espritz_features(
    sequence_rows: Iterable[dict[str, object]],
    *,
    cache_path: str | Path,
    model: str = "D",
    sw: int = 0,
    chunk_size: int = 250,
    workers: int = 1,
) -> list[ESpritzFeature]:
    cache = ESpritzFeatureCache(cache_path)
    cached_features = cache._read_all()
    rows = list(sequence_rows)
    features: list[ESpritzFeature] = []
    pending_rows: list[dict[str, object]] = []
    for row in rows:
        accession = str(row["UniprotEntry"]).strip()
        cached = cached_features.get(accession)
        if cached is not None and cached.model == model and cached.sw == sw:
            features.append(cached)
            continue
        pending_rows.append(row)

    chunks = [
        pending_rows[start : start + chunk_size]
        for start in range(0, len(pending_rows), chunk_size)
    ]
    if workers <= 1:
        for chunk in chunks:
            for feature in compute_espritz_features(chunk, model=model, sw=sw):
                cache.append_new(feature)
                cached_features[feature.accession] = feature
                features.append(feature)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(compute_espritz_features, chunk, model=model, sw=sw)
                for chunk in chunks
            ]
            for future in as_completed(futures):
                for feature in future.result():
                    cache.append_new(feature)
                    cached_features[feature.accession] = feature
                    features.append(feature)

    by_accession = {feature.accession: feature for feature in features}
    return [
        by_accession[str(row["UniprotEntry"]).strip()]
        for row in rows
        if str(row["UniprotEntry"]).strip() in by_accession
    ]


def compute_espritz_features(
    sequence_rows: Iterable[dict[str, object]],
    *,
    model: str = "D",
    sw: int = 0,
) -> list[ESpritzFeature]:
    if not ESPRITZ_RUNNER.exists():
        raise FileNotFoundError(f"ESpritz runner not found: {ESPRITZ_RUNNER}")
    rows = list(sequence_rows)
    with tempfile.TemporaryDirectory(prefix="phasepred_espritz_") as tmp:
        workdir = Path(tmp)
        accessions: list[str] = []
        for row in rows:
            accession = str(row["UniprotEntry"]).strip()
            sequence = "".join(str(row["sequence"]).split()).upper()
            if not accession or not sequence:
                continue
            accessions.append(accession)
            (workdir / f"{_safe_stem(accession)}.fasta").write_text(
                f">{accession}\n{sequence}\n",
                encoding="utf-8",
            )
        if not accessions:
            return []
        subprocess.run(
            [str(ESPRITZ_RUNNER), str(workdir), model, str(sw)],
            check=True,
            cwd=ESPRITZ_RUNNER.parents[2],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        features = []
        for accession in accessions:
            output = workdir / f"{_safe_stem(accession)}.espritz"
            if output.exists():
                features.append(parse_espritz_output(output, accession, model=model, sw=sw))
        return features


def parse_espritz_output(
    path: str | Path,
    accession: str,
    *,
    model: str = "D",
    sw: int = 0,
) -> ESpritzFeature:
    states: list[str] = []
    scores: list[float] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.strip().split()
        if len(parts) != 2 or parts[0] not in {"D", "O"}:
            continue
        try:
            score = float(parts[1])
        except ValueError:
            continue
        states.append(parts[0])
        scores.append(score)
    if not scores:
        raise ValueError(f"No ESpritz residue scores parsed from {path}")
    return ESpritzFeature(
        accession=accession,
        espritz_idr=sum(state == "D" for state in states) / len(states),
        mean_score=sum(scores) / len(scores),
        max_score=max(scores),
        disorder_threshold=0.5072 if model == "D" and sw == 0 else float("nan"),
        residues=len(scores),
        model=model,
        sw=sw,
    )


def _safe_stem(accession: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in accession)
