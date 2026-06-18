from __future__ import annotations

import json
import sys
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

IUPRED_DISORDER_THRESHOLD = 0.5


def _iupred_tool_dir() -> Path:
    """Resolve the IUPred3 install dir via tool_paths (env / vendored / PATH)."""
    from phasepred.tool_paths import find_iupred3_dir
    return find_iupred3_dir()


@dataclass(frozen=True)
class IUPredFeature:
    accession: str
    iupred_idr: float
    mean_score: float
    max_score: float
    disorder_threshold: float
    residues: int
    mode: str
    smoothing: str

    def to_json(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_json(cls, payload: dict[str, object]) -> IUPredFeature:
        return cls(
            accession=str(payload["accession"]),
            iupred_idr=float(payload["iupred_idr"]),
            mean_score=float(payload["mean_score"]),
            max_score=float(payload["max_score"]),
            disorder_threshold=float(payload["disorder_threshold"]),
            residues=int(payload["residues"]),
            mode=str(payload["mode"]),
            smoothing=str(payload["smoothing"]),
        )


class IUPredFeatureCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def get(self, accession: str) -> IUPredFeature | None:
        return self._read_all().get(accession.strip())

    def put(self, feature: IUPredFeature) -> None:
        features = self._read_all()
        if feature.accession not in features:
            self.append_new(feature)
            return
        features[feature.accession] = feature
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("w", encoding="utf-8") as handle:
            for item in features.values():
                handle.write(json.dumps(item.to_json(), sort_keys=True) + "\n")

    def append_new(self, feature: IUPredFeature) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(feature.to_json(), sort_keys=True) + "\n")

    def _read_all(self) -> dict[str, IUPredFeature]:
        if not self.path.exists():
            return {}
        features: dict[str, IUPredFeature] = {}
        with self.path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                feature = IUPredFeature.from_json(json.loads(stripped))
                features[feature.accession] = feature
        return features


def compute_iupred_feature(
    accession: str,
    sequence: str,
    *,
    mode: str = "long",
    smoothing: str = "medium",
    threshold: float = IUPRED_DISORDER_THRESHOLD,
) -> IUPredFeature:
    scores = compute_iupred_scores(sequence, mode=mode, smoothing=smoothing)
    if not scores:
        raise ValueError(f"IUPred3 returned no scores for {accession}")
    return IUPredFeature(
        accession=accession,
        iupred_idr=sum(score >= threshold for score in scores) / len(scores),
        mean_score=sum(scores) / len(scores),
        max_score=max(scores),
        disorder_threshold=threshold,
        residues=len(scores),
        mode=mode,
        smoothing=smoothing,
    )


def compute_iupred_scores(
    sequence: str,
    *,
    mode: str = "long",
    smoothing: str = "medium",
) -> list[float]:
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise ValueError("Cannot compute IUPred3 scores for an empty sequence")
    module = _load_iupred_module()
    scores, _ = module.iupred(normalized, mode, smoothing=smoothing)
    return [float(score) for score in scores]


def resolve_iupred_features(
    sequence_rows: Iterable[dict[str, object]],
    *,
    cache_path: str | Path,
    mode: str = "long",
    smoothing: str = "medium",
    threshold: float = IUPRED_DISORDER_THRESHOLD,
    workers: int = 1,
) -> list[IUPredFeature]:
    cache = IUPredFeatureCache(cache_path)
    cached_features = cache._read_all()
    features: list[IUPredFeature] = []
    pending_rows = []
    for row in sequence_rows:
        accession = str(row["UniprotEntry"]).strip()
        cached = cached_features.get(accession)
        if cached is not None:
            features.append(cached)
            continue
        pending_rows.append((accession, str(row["sequence"])))
    if workers <= 1:
        for accession, sequence in pending_rows:
            feature = compute_iupred_feature(
                accession,
                sequence,
                mode=mode,
                smoothing=smoothing,
                threshold=threshold,
            )
            cache.put(feature)
            cached_features[accession] = feature
            features.append(feature)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    _compute_iupred_feature_worker,
                    accession,
                    sequence,
                    mode,
                    smoothing,
                    threshold,
                ): accession
                for accession, sequence in pending_rows
            }
            for future in as_completed(futures):
                feature = future.result()
                cache.append_new(feature)
                cached_features[feature.accession] = feature
                features.append(feature)
    return features


def _compute_iupred_feature_worker(
    accession: str,
    sequence: str,
    mode: str,
    smoothing: str,
    threshold: float,
) -> IUPredFeature:
    return compute_iupred_feature(
        accession,
        sequence,
        mode=mode,
        smoothing=smoothing,
        threshold=threshold,
    )


def _load_iupred_module():
    tool_dir = _iupred_tool_dir()
    tool_path = str(tool_dir)
    if tool_path not in sys.path:
        sys.path.insert(0, tool_path)
    import iupred3_lib

    return iupred3_lib
