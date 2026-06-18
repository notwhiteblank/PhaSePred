"""PhaSePred predictor: compute features and score proteins from FASTA."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd


class PhaSePredMissingDataWarning(UserWarning):
    """Emitted when an optional data file (DeepPhase, PhosphoSitePlus) is
    missing or does not cover specific accessions, so the corresponding
    feature column will fall back to median imputation in the pipeline."""


# Repo root relative to this module
REPO_ROOT = Path(__file__).resolve().parents[2]

# Feature columns
from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS

# Tool path resolver (handles ENV → vendored → PATH cascade)
from phasepred.tool_paths import (
    PhaSePredToolNotFound,
    assert_predict_requirements,
    find_espritz_wrapper,
    find_plaac_wrapper,
    find_pscore_dir,
    find_pscore_script,
    find_seg,
)

# Cache paths
DEEPHASE_TABLE = REPO_ROOT / "data" / "raw" / "external" / "deepphase" / "extracted" / "tableS3.xlsx"
PHOSPHOSITE_PATH = REPO_ROOT / "data" / "raw" / "external" / "phosphositeplus" / "Phosphorylation_site_dataset.gz"

# catGRANULE v1 Z-normalization constants — bundled with the package because
# they are model calibration parameters (mean, std of training-set raw scores)
# rather than user data. Loaded from src/phasepred/data/catgranule_v1_norm.json.
CATGRANULE_NORM = Path(__file__).parent / "data" / "catgranule_v1_norm.json"

ESPRITZ_CACHE = REPO_ROOT / "data" / "interim" / "espritz_idr_cache.jsonl"

# Default product directory. Product B (extended dataset) is the default
# because its held-out AUC is higher; Product A reproduces the paper split.
DEFAULT_PRODUCT = "B"
PRODUCT_MODEL_DIRS = {
    "A": REPO_ROOT / "products" / "A_paper_split_recomputed" / "models",
    "B": REPO_ROOT / "products" / "B_extended_dataset" / "models",
}
DEFAULT_MODELS_DIR = PRODUCT_MODEL_DIRS[DEFAULT_PRODUCT]


def _safe_stem(accession: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in accession)


# ---------------------------------------------------------------------------
# Feature computation
# ---------------------------------------------------------------------------


def compute_all_features(
    records: list[dict[str, str]],
    *,
    include_human: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """Compute all PhaSePred features for a list of FASTA records.

    Each record is {"accession": str, "sequence": str}.
    Returns a DataFrame with feature columns ready for prediction.
    """
    t_start = time.time()

    # --- Native features (pure Python, fast) ---
    from phasepred.features import compute_native_features as _native

    native_rows = []
    for r in records:
        seq = str(r["sequence"]).upper()
        try:
            feats = _native(seq)
            native_rows.append(
                {
                    "UniprotEntry": r["accession"],
                    "length": feats["length"],
                    "Hydropathy": feats["Hydropathy"],
                    "FCR": feats["FCR"],
                }
            )
        except Exception:
            native_rows.append(
                {
                    "UniprotEntry": r["accession"],
                    "length": float("nan"),
                    "Hydropathy": float("nan"),
                    "FCR": float("nan"),
                }
            )

    feature_df = pd.DataFrame(native_rows)
    if verbose:
        print(f"  Native features: {len(feature_df)} in {time.time() - t_start:.1f}s")

    # --- catGRANULE v1 (pure Python) ---
    t0 = time.time()
    from phasepred.catgranule_v1 import CatGranuleNormalization, normalize_score, raw_score

    cat_scores: dict[str, float] = {}
    if CATGRANULE_NORM.exists():
        with open(CATGRANULE_NORM) as fh:
            norm_data = json.load(fh)
        norm = CatGranuleNormalization(
            mean=float(norm_data["mean"]),
            standard_deviation=float(norm_data["standard_deviation"]),
        )
        for r in records:
            try:
                raw = raw_score(str(r["sequence"]))
                cat_scores[r["accession"]] = normalize_score(raw, norm)
            except Exception:
                cat_scores[r["accession"]] = float("nan")
    feature_df["catGRANULE"] = feature_df["UniprotEntry"].map(cat_scores)
    if verbose:
        n = feature_df["catGRANULE"].notna().sum()
        print(f"  catGRANULE v1: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- ESpritz IDR ---
    t0 = time.time()
    espritz_scores = _run_espritz(records)
    feature_df["IDR"] = feature_df["UniprotEntry"].map(espritz_scores)
    if verbose:
        n = feature_df["IDR"].notna().sum()
        print(f"  ESpritz IDR: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- SEG LCR ---
    t0 = time.time()
    seg_scores = _run_seg(records)
    feature_df["LCR"] = feature_df["UniprotEntry"].map(seg_scores)
    if verbose:
        n = feature_df["LCR"].notna().sum()
        print(f"  SEG LCR: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- PScore (parallel) ---
    t0 = time.time()
    pscore_scores = _run_pscore_parallel(records)
    feature_df["PScore"] = feature_df["UniprotEntry"].map(pscore_scores)
    if verbose:
        n = feature_df["PScore"].notna().sum()
        print(f"  PScore: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- PLAAC (parallel) ---
    t0 = time.time()
    plaac_scores = _run_plaac_parallel(records)
    feature_df["PLAAC"] = feature_df["UniprotEntry"].map(plaac_scores)
    if verbose:
        n = feature_df["PLAAC"].notna().sum()
        print(f"  PLAAC: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- DeepCoil (on-demand subprocess into isolated conda env) ---
    t0 = time.time()
    dc_scores = _run_deepcoil(records)
    feature_df["DeepCoil"] = feature_df["UniprotEntry"].map(dc_scores)
    if verbose:
        n = feature_df["DeepCoil"].notna().sum()
        print(f"  DeepCoil: {n}/{len(feature_df)} in {time.time() - t0:.1f}s")

    # --- Human-specific features ---
    if include_human:
        t0 = time.time()
        phos_scores = _compute_phos_freq(records)
        feature_df["Phos freq"] = feature_df["UniprotEntry"].map(phos_scores)
        dp_scores = _lookup_deepphase([r["accession"] for r in records])
        feature_df["DeepPhase"] = feature_df["UniprotEntry"].map(dp_scores)
        if verbose:
            n_p = feature_df["Phos freq"].notna().sum()
            n_d = feature_df["DeepPhase"].notna().sum()
            print(f"  Phos freq + DeepPhase: {n_p}, {n_d} in {time.time() - t0:.1f}s")

    if verbose:
        total_t = time.time() - t_start
        print(f"  Total feature computation: {total_t:.0f}s")

    return feature_df


# ---------------------------------------------------------------------------
# ESpritz runner
# ---------------------------------------------------------------------------


def _run_espritz(records: list[dict[str, str]]) -> dict[str, float]:
    """Run ESpritz DisProt 5% FPR on all records (batch in one directory)."""
    workdir = Path(tempfile.mkdtemp(prefix="espritz_pred_"))
    try:
        for r in records:
            fasta = workdir / f"{_safe_stem(r['accession'])}.fasta"
            fasta.write_text(
                f">{r['accession']}\n{str(r['sequence']).upper()}\n", encoding="utf-8"
            )
        result = subprocess.run(
            [str(find_espritz_wrapper()), str(workdir), "D", "0"],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            return {}

        scores: dict[str, float] = {}
        for r in records:
            out = workdir / f"{_safe_stem(r['accession'])}.espritz"
            if out.exists():
                text = out.read_text(encoding="utf-8")
                lines = [l.strip() for l in text.split("\n") if l.strip() and not l.startswith("#")]
                if lines:
                    states = "".join(lines)
                    scores[r["accession"]] = sum(s == "D" for s in states) / len(states)
        return scores
    finally:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# SEG runner
# ---------------------------------------------------------------------------


def _run_seg(records: list[dict[str, str]]) -> dict[str, float]:
    """Run SEG on all records."""
    try:
        seg_bin = find_seg()
    except PhaSePredToolNotFound:
        return {}

    fasta_text = "".join(
        f">{r['accession']}\n{str(r['sequence']).upper()}\n" for r in records
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, prefix="seg_pred_") as fh:
        fh.write(fasta_text)
        fasta_path = Path(fh.name)

    try:
        result = subprocess.run(
            [str(seg_bin), str(fasta_path), "12", "2.2", "2.5", "-x"],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            return {}

        masked: dict[str, str] = {}
        current_id = ""
        current_seq: list[str] = []
        for line in result.stdout.split("\n"):
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    masked[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_id:
            masked[current_id] = "".join(current_seq)

        scores: dict[str, float] = {}
        for r in records:
            seq = masked.get(r["accession"], "")
            if seq:
                lc = sum(1 for c in seq if c == "x" or c.islower())
                scores[r["accession"]] = lc / len(seq)
        return scores
    finally:
        fasta_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# PScore (parallel)
# ---------------------------------------------------------------------------


def _pscore_chunk(args: tuple[int, list[dict[str, str]]]) -> tuple[int, str]:
    chunk_idx, recs = args
    fasta_text = "".join(
        f">{r['accession']}\n{str(r['sequence']).upper()}\n" for r in recs
    )
    with tempfile.TemporaryDirectory(prefix=f"pscore_pred_{chunk_idx}_") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "input.fasta"
        output = workdir / "pscore.tsv"
        fasta.write_text(fasta_text, encoding="utf-8")
        subprocess.run(
            ["uv", "run", "python", str(find_pscore_script()),
             str(fasta), "-output", str(output), "-overwrite", "-mute"],
            cwd=str(find_pscore_dir()),
            check=True, capture_output=True, text=True, timeout=300,
        )
        return (chunk_idx, output.read_text(encoding="utf-8"))


def _run_pscore_parallel(records: list[dict[str, str]], n_workers: int = 8) -> dict[str, float]:
    """Run PScore in parallel chunks."""
    try:
        find_pscore_script()
    except PhaSePredToolNotFound:
        return {}

    chunk_size = max(1, len(records) // n_workers)
    chunks = []
    for i in range(0, len(records), chunk_size):
        chunks.append((len(chunks), records[i : i + chunk_size]))

    scores: dict[str, float] = {}
    with Pool(processes=min(n_workers, len(chunks))) as pool:
        for _, output in pool.imap_unordered(_pscore_chunk, chunks):
            for line in output.strip().split("\n"):
                if not line.startswith("PScore:"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        value = float(parts[1])
                        acc = parts[-1].lstrip(">")
                        scores[acc] = value
                    except (ValueError, IndexError):
                        continue
    return scores


# ---------------------------------------------------------------------------
# PLAAC (parallel)
# ---------------------------------------------------------------------------


def _plaac_chunk(args: tuple[int, list[dict[str, str]]]) -> tuple[int, str]:
    chunk_idx, recs = args
    fasta_text = "".join(
        f">{r['accession']}\n{str(r['sequence']).upper()}\n" for r in recs
    )
    with tempfile.TemporaryDirectory(prefix=f"plaac_pred_{chunk_idx}_") as tmp:
        fasta = Path(tmp) / "input.fasta"
        fasta.write_text(fasta_text, encoding="utf-8")
        result = subprocess.run(
            [str(find_plaac_wrapper()), "-i", str(fasta)],
            cwd=str(REPO_ROOT), check=True, capture_output=True, text=True, timeout=300,
        )
        return (chunk_idx, result.stdout)


def _run_plaac_parallel(records: list[dict[str, str]], n_workers: int = 4) -> dict[str, float]:
    """Run PLAAC in parallel chunks."""
    try:
        find_plaac_wrapper()
    except PhaSePredToolNotFound:
        return {}

    chunk_size = max(1, len(records) // n_workers)
    chunks = []
    for i in range(0, len(records), chunk_size):
        chunks.append((len(chunks), records[i : i + chunk_size]))

    scores: dict[str, float] = {}
    with Pool(processes=min(n_workers, len(chunks))) as pool:
        for _, output in pool.imap_unordered(_plaac_chunk, chunks):
            for line in output.strip().split("\n"):
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0] != "ACCESSION":
                    try:
                        scores[str(parts[0])] = float(parts[2])
                    except (ValueError, IndexError):
                        continue
    return scores


# ---------------------------------------------------------------------------
# Lookup / runner helpers
# ---------------------------------------------------------------------------


DEEPCOIL_WRAPPER = REPO_ROOT / "tools" / "wrappers" / "run_deepcoil.sh"


def _safe_id(accession: str) -> str:
    """DeepCoil sanitizes output filenames; mirror its rule."""
    return "".join(ch for ch in str(accession) if ch.isalnum() or ch == "_")


def _summarize_deepcoil_output(out_path: Path) -> float:
    """Parse a DeepCoil per-residue .out file → mean raw_cc."""
    frame = pd.read_csv(out_path, sep="\t")
    if frame.empty:
        raise ValueError("empty DeepCoil output")
    return float(frame["raw_cc"].mean())


def _run_deepcoil(records: list[dict[str, str]]) -> dict[str, float]:
    """Run DeepCoil on demand for the user's input sequences.

    Aggregates per-residue raw_cc into a single per-protein mean — this is
    the DeepCoil feature definition used by both products (matches the
    pre-computed reference table the project used during training).

    Returns {} on any failure (wrapper not installed, DeepCoil env broken,
    subprocess timeout). The caller already tolerates NaN columns via the
    pipeline's SimpleImputer; this function does not warn (DeepCoil being
    unavailable is a tool-installation issue surfaced by check-tools).
    """
    if not DEEPCOIL_WRAPPER.is_file():
        return {}
    if not records:
        return {}

    scores: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="deepcoil_predict_") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "input.fasta"
        out_dir = workdir / "out"
        out_dir.mkdir()

        lines: list[str] = []
        id_map: dict[str, str] = {}  # safe_id → original accession
        for r in records:
            acc = str(r["accession"])
            safe = _safe_id(acc)
            id_map[safe] = acc
            seq = re.sub(r"\s+", "", str(r["sequence"]).upper())
            lines.append(f">{safe}")
            lines.extend(seq[i : i + 80] for i in range(0, len(seq), 80))
        fasta.write_text("\n".join(lines) + "\n", encoding="utf-8")

        try:
            subprocess.run(
                [
                    str(DEEPCOIL_WRAPPER),
                    "-i", str(fasta),
                    "-out_path", str(out_dir),
                    "-n_cpu", "4",
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=max(300, 30 * len(records)),
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return {}

        for safe, orig in id_map.items():
            out_file = out_dir / f"{safe}.out"
            if not out_file.exists():
                continue
            try:
                scores[orig] = _summarize_deepcoil_output(out_file)
            except Exception:
                continue
    return scores


def _lookup_deepphase(accessions: list[str]) -> dict[str, float]:
    """Look up DeepPhase scores from upstream table S3.

    DeepPhase has no public inference tool — only the supplement table for
    ~6000 reference human proteins. Proteins not in the table can't be
    scored; we warn the user so they see which accessions are imputed.
    """
    if not DEEPHASE_TABLE.exists():
        warnings.warn(
            (
                "DeepPhase table missing at "
                f"{DEEPHASE_TABLE.relative_to(REPO_ROOT) if str(DEEPHASE_TABLE).startswith(str(REPO_ROOT)) else DEEPHASE_TABLE}. "
                "Human-mode predictions (hSaPS / hPdPS) will use median "
                "imputation for the DeepPhase column. See the README "
                "section 'Optional human-feature data' for the download."
            ),
            PhaSePredMissingDataWarning,
            stacklevel=2,
        )
        return {}
    dp = pd.read_excel(str(DEEPHASE_TABLE), sheet_name="tableS3")
    dp_map: dict[str, float] = {}
    for _, row in dp.iterrows():
        swiss = str(row.get("Swiss-Prot ID", "")).strip()
        if swiss and swiss != "nan":
            try:
                dp_map[swiss] = float(row["DeepPhase_score"])
            except (ValueError, KeyError):
                continue
    hits = {a: dp_map[a] for a in accessions if a in dp_map}
    missing = [a for a in accessions if a not in dp_map]
    if missing:
        preview = ", ".join(missing[:5]) + ("…" if len(missing) > 5 else "")
        warnings.warn(
            (
                f"DeepPhase has no entry for {len(missing)}/{len(accessions)} "
                f"requested protein(s): {preview}. Median imputation will be "
                "used, so hSaPS/hPdPS scores for these proteins should be "
                "treated as approximate."
            ),
            PhaSePredMissingDataWarning,
            stacklevel=2,
        )
    return hits


def _compute_phos_freq(records: list[dict[str, str]]) -> dict[str, float]:
    """Compute phosphorylation site frequency from PhosphoSitePlus.

    Cannot be re-computed without the PhosphoSitePlus dataset, which the
    upstream license forbids redistributing. We warn loudly when the file
    is missing or specific proteins aren't covered.
    """
    if not PHOSPHOSITE_PATH.exists():
        warnings.warn(
            (
                "PhosphoSitePlus dataset missing at "
                f"{PHOSPHOSITE_PATH.relative_to(REPO_ROOT) if str(PHOSPHOSITE_PATH).startswith(str(REPO_ROOT)) else PHOSPHOSITE_PATH}. "
                "Human-mode predictions (hSaPS / hPdPS) will use median "
                "imputation for the 'Phos freq' column. See the README "
                "section 'Optional human-feature data' for the registration "
                "and download steps."
            ),
            PhaSePredMissingDataWarning,
            stacklevel=2,
        )
        return {}
    from phasepred.legacy_features import compute_phosphosite_frequencies

    seq_rows = [{"UniprotEntry": r["accession"], "sequence": r["sequence"]} for r in records]
    features = compute_phosphosite_frequencies(seq_rows, PHOSPHOSITE_PATH)
    hits = {f.accession: f.value for f in features}
    requested = {r["accession"] for r in records}
    # compute_phosphosite_frequencies only emits accessions that appeared in
    # the PhosphoSitePlus table; everything else is silently dropped.
    # A row with value 0.0 is fine (no phospho sites for that protein);
    # absence means the protein isn't in the database.
    missing = sorted(requested - set(hits.keys()))
    if missing:
        preview = ", ".join(missing[:5]) + ("…" if len(missing) > 5 else "")
        warnings.warn(
            (
                f"PhosphoSitePlus has no entry for {len(missing)}/{len(requested)} "
                f"requested protein(s): {preview}. Median imputation will be "
                "used for 'Phos freq', so hSaPS/hPdPS scores for these "
                "proteins should be treated as approximate."
            ),
            PhaSePredMissingDataWarning,
            stacklevel=2,
        )
    return hits


# ---------------------------------------------------------------------------
# Model loading and prediction
# ---------------------------------------------------------------------------


def load_ensemble(model_dir: Path) -> list:
    """Load a 10-model XGBoost ensemble from a directory of joblib files."""
    import joblib

    model_files = sorted(model_dir.glob("8f_model_*.joblib"))
    if not model_files:
        raise FileNotFoundError(f"No 8f_model_*.joblib files found in {model_dir}")
    return [joblib.load(str(f)) for f in model_files]


def load_single_model(model_path: Path):
    """Load a single model artifact (joblib file)."""
    import joblib

    if not model_path.exists():
        raise FileNotFoundError(f"Model not found: {model_path}")
    return joblib.load(str(model_path))


def predict_scores(
    models: list,
    features: pd.DataFrame,
    feature_columns: list[str],
) -> pd.DataFrame:
    """Predict scores using a loaded ensemble (or single model).

    If the model was trained on a subset of columns (e.g. paper baseline
    with only 3 features), use only those columns from the feature table.

    Returns DataFrame with columns: UniprotEntry, score
    """
    # Check if the model artifact specifies its own feature columns
    model_columns = None
    if isinstance(models, list) and models:
        first = models[0]
        if hasattr(first, "feature_columns"):
            model_columns = list(first.feature_columns)
    elif hasattr(models, "feature_columns"):
        model_columns = list(models.feature_columns)

    effective_columns = model_columns or feature_columns
    avail_cols = [c for c in effective_columns if c in features.columns]
    missing = set(effective_columns) - set(avail_cols)
    if missing:
        for col in missing:
            features[col] = float("nan")

    X = features[avail_cols]

    # models may be a single object or a list
    model_list = models if isinstance(models, list) else [models]

    if len(model_list) > 1:
        scores = np.zeros(len(X))
        for model in model_list:
            scores += _model_predict_proba(model, X)[:, 1]
        scores /= len(model_list)
    else:
        scores = _model_predict_proba(model_list[0], X)[:, 1]

    result = features[["UniprotEntry"]].copy()
    result["score"] = scores
    return result


def _model_predict_proba(model, X: pd.DataFrame) -> np.ndarray:
    """Call predict_proba, handling both sklearn Pipeline and ModelArtifact."""
    if hasattr(model, "estimator"):
        # ModelArtifact wrapper
        return model.estimator.predict_proba(X)
    return model.predict_proba(X)
