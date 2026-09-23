"""External feature-tool runners.

Module boundary (S2): ``features`` holds pure computation, while this
module owns the *four-step* lifecycle for every external feature tool —
**locate** (resolve the tool-package entry), **execute** (subprocess or
lookup), **parse** (translate raw output into per-accession values), and
**degrade** (return ``{}`` so the caller stamps NaN).

Since S3 the *locate* step goes through the manifest registry in
``phasepred.tool_paths``: each tool is a self-contained package under
``tools/<TOOL>/`` (``run`` entrypoint + ``manifest.toml`` contract) and the
registry resolves ``PHASEPRED_<TOOL>_DIR`` → ``tools/<TOOL>/`` → PATH.
The package ``run`` scripts dispatch exactly the same underlying commands
as the pre-S3 wrappers, so output values are unchanged.

The universal public surface is ``run_<tool>(records) -> dict[str, float]``
with ``records`` a list of ``{"accession": str, "sequence": str}`` dicts.
A missing/broken tool degrades to ``{}`` — never raises at inference time
(path resolution problems still raise ``PhaSePredToolNotFound`` so
``predict``'s preflight can fail fast).

Per-feature v2022 definitions (``docs/FEATURES.md``):
- DeepCoil: binarized ``1.0 if max(raw_cc) >= 0.82 else 0.0``.
- ESpritz: DisProt model ``D``, ``sw 0`` (5% FPR), IDR = fraction of ``D``
  state residues parsed from the two-column residue lines only.
"""

from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
import warnings
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from phasepred.data import bundled_path, data_path, optional_data_path
from phasepred.tool_paths import (
    contract_env,
    find_tool_entry,
    find_tool_entry_or_none,
)


class PhaSePredMissingDataWarning(UserWarning):
    """Emitted when optional data (DeepPhase/PhosphoSitePlus) is missing.

    The lookup runners are the only emitters; ``predictor`` re-exports it so
    ``predict`` can register it with the ``warnings`` machinery.
    """


def deepphase_table_path() -> Path:
    """Resolve the DeepPhase score table lazily (at call time).

    Preference: the package-bundled derived TSV shipped in the wheel
    (``data/deepphase/deepphase_scores.tsv``) > the optional external xlsx
    (``optional_data_path``). ``lookup_deepphase`` reads whichever this
    returns and degrades to ``{}`` when the resolved path does not exist.
    """
    bundled = bundled_path("data", "deepphase", "deepphase_scores.tsv")
    if bundled.exists():
        return bundled
    return optional_data_path("deepphase", "extracted", "tableS3.xlsx")


def phosphosite_path() -> Path:
    """Resolve the PhosphoSitePlus dataset lazily (at call time)."""
    return optional_data_path("phosphositeplus", "Phosphorylation_site_dataset.gz")


def phosphosite_version_marker(path: Path | None = None) -> str | None:
    """Return the local PhosphoSitePlus version marker, or ``None``.

    The marker is the first non-empty line of the gunzipped dataset (e.g.
    ``"042026"``; see ``tools/PhosphoSitePlus/manifest.toml``). Ruling 2: any
    comparison against a frozen ``Phos freq`` value must first read this
    marker -- a match means the frozen value can be compared exactly, a
    mismatch must SKIP rather than fail, because PhosphoSitePlus is a rolling
    database (upstream had already rolled from ``042026`` to ``071726`` during
    E4).
    """
    dataset = path if path is not None else phosphosite_path()
    if not dataset.is_file():
        return None
    try:
        with gzip.open(dataset, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                marker = line.strip()
                if marker:
                    return marker
    except OSError:
        return None
    return None


ESPRITZ_STATE_RE = re.compile(r"^(?P<state>[DO])\s+(?P<score>\S+)$")

# ---------------------------------------------------------------------------
# ESpritz (IDR)
# ---------------------------------------------------------------------------


def _safe_stem(accession: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in accession)


def run_espritz(
    records: list[dict[str, str]],
    *,
    model: str = "D",
    sw: int = 0,
) -> dict[str, float]:
    """Run ESpritz DisProt 5% FPR on all records (batch in one directory).

    Returns the per-accession IDR fraction (``count(D) / L``) computed from
    the two-column residue lines only — banner/license lines are ignored.
    """
    wrapper = find_tool_entry("espritz")
    workdir = Path(tempfile.mkdtemp(prefix="espritz_pred_"))
    try:
        for r in records:
            fasta = workdir / f"{_safe_stem(r['accession'])}.fasta"
            fasta.write_text(
                f">{r['accession']}\n{str(r['sequence']).upper()}\n", encoding="utf-8"
            )
        result = subprocess.run(
            [str(wrapper), str(workdir), model, str(sw)],
            capture_output=True, text=True, timeout=300,
            env=contract_env(),
        )
        if result.returncode != 0:
            return {}

        scores: dict[str, float] = {}
        for r in records:
            out = workdir / f"{_safe_stem(r['accession'])}.espritz"
            if not out.exists():
                continue
            states = _parse_espritz_states(out)
            if states:
                scores[r["accession"]] = sum(s == "D" for s in states) / len(states)
        return scores
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _parse_espritz_states(out_path: Path) -> list[str]:
    """Extract ordered residue states from a ``.espritz`` output file.

    Only two-column residue lines whose first token is D or O are counted;
    license/banner/sequence lines are skipped (v2022 IDR parsing rule).
    """
    states: list[str] = []
    for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ESPRITZ_STATE_RE.match(line.strip())
        if match is not None:
            states.append(match.group("state"))
    return states


# ---------------------------------------------------------------------------
# SEG (LCR)
# ---------------------------------------------------------------------------


def run_seg(records: list[dict[str, str]]) -> dict[str, float]:
    """Run SEG with the default parameters (``12 2.2 2.5 -x``)."""
    seg_bin = find_tool_entry("seg")

    fasta_text = "".join(
        f">{r['accession']}\n{str(r['sequence']).upper()}\n" for r in records
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="seg_pred_"
    ) as fh:
        fh.write(fasta_text)
        fasta_path = Path(fh.name)

    try:
        result = subprocess.run(
            [str(seg_bin), str(fasta_path), "12", "2.2", "2.5", "-x"],
            capture_output=True, text=True, timeout=120,
            env=contract_env(),
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
# PScore
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
            [str(find_tool_entry("pscore")),
             str(fasta), "-output", str(output), "-overwrite", "-mute"],
            check=True, capture_output=True, text=True, timeout=300,
            env=contract_env(),
        )
        return (chunk_idx, output.read_text(encoding="utf-8"))


def run_pscore(
    records: list[dict[str, str]],
    *,
    n_workers: int = 8,
) -> dict[str, float]:
    """Run PScore in parallel chunks; single values from ``PScore:`` lines."""
    find_tool_entry("pscore")

    chunk_size = max(1, len(records) // n_workers)
    chunks: list[tuple[int, list[dict[str, str]]]] = []
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
# PLAAC (NLLR)
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
            [str(find_tool_entry("plaac")), "-i", str(fasta)],
            cwd=str(data_path()),
            check=True, capture_output=True, text=True, timeout=300,
            env=contract_env(),
        )
        return (chunk_idx, result.stdout)


def run_plaac(
    records: list[dict[str, str]],
    *,
    n_workers: int = 4,
) -> dict[str, float]:
    """Run PLAAC in parallel chunks; take the ``NLLR`` column value."""
    find_tool_entry("plaac")

    chunk_size = max(1, len(records) // n_workers)
    chunks: list[tuple[int, list[dict[str, str]]]] = []
    for i in range(0, len(records), chunk_size):
        chunks.append((len(chunks), records[i : i + chunk_size]))

    scores: dict[str, float] = {}
    with Pool(processes=min(n_workers, len(chunks))) as pool:
        for _, output in pool.imap_unordered(_plaac_chunk, chunks):
            header: list[str] | None = None
            for line in output.strip().split("\n"):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    if parts and parts[0] == "SEQid":
                        header = parts
                    continue
                if len(parts) != len(header):
                    continue
                record = dict(zip(header, parts, strict=False))
                try:
                    scores[str(record["SEQid"])] = float(record["NLLR"])
                except (KeyError, ValueError):
                    continue
    return scores


# ---------------------------------------------------------------------------
# DeepCoil (v2022 binary definition)
# ---------------------------------------------------------------------------


def _safe_id(accession: str) -> str:
    """DeepCoil sanitizes output filenames; mirror its rule."""
    return "".join(ch for ch in str(accession) if ch.isalnum() or ch == "_")


DEEPCOIL_THRESHOLD = 0.82


def _binarize_deepcoil_frame(frame: pd.DataFrame) -> float:
    """v2022 DeepCoil feature: ``1.0 if max(raw_cc) >= 0.82 else 0.0``."""
    if frame.empty:
        raise ValueError("empty DeepCoil output")
    return float(1.0 if frame["raw_cc"].max() >= DEEPCOIL_THRESHOLD else 0.0)


def run_deepcoil(records: list[dict[str, str]]) -> dict[str, float]:
    """Run DeepCoil on demand for the input sequences.

    Aggregates the per-residue ``raw_cc`` array into the v2022 *binary*
    coiled-coil indicator (threshold 0.82). Returns ``{}`` on any failure
    (wrapper not installed, env broken, subprocess timeout) so the caller
    degrades the column to NaN.
    """
    if not records:
        return {}

    runner = find_tool_entry_or_none("deepcoil")
    if runner is None:
        return {}  # contract package missing → degrade to NaN

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
                    str(runner),
                    "-i", str(fasta),
                    "-out_path", str(out_dir),
                    "-n_cpu", "4",
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=max(300, 30 * len(records)),
                env=contract_env(),
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return {}

        for safe, orig in id_map.items():
            out_file = out_dir / f"{safe}.out"
            if not out_file.exists():
                continue
            try:
                frame = pd.read_csv(out_file, sep="\t")
                scores[orig] = _binarize_deepcoil_frame(frame)
            except Exception:
                continue
    return scores


# ---------------------------------------------------------------------------
# Residue-level extractors (tool-features subcommand)
#
# These mirror the *per-residue* arrays of the web ``tools.<x>`` payloads for
# the tools whose local contracts expose them. Single values stay on the
# ``run_<tool>`` surface; here we return dict[accession, dict] with the raw
# residue lists so the CLI can derive both the summary and the JSON.
# ---------------------------------------------------------------------------


def _parse_espritz_residues(out_path: Path) -> tuple[list[str], list[float]]:
    states: list[str] = []
    scores: list[float] = []
    for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ESPRITZ_STATE_RE.match(line.strip())
        if match is not None:
            states.append(match.group("state"))
            try:
                scores.append(float(match.group("score")))
            except ValueError:
                scores.append(float("nan"))
    return states, scores


def espritz_residues(
    records: list[dict[str, str]],
    *,
    model: str = "D",
    sw: int = 0,
) -> dict[str, dict[str, object]]:
    """Per-residue ESpritz ``label`` (D/O) + ``residue`` (score) arrays."""
    wrapper = find_tool_entry("espritz")
    workdir = Path(tempfile.mkdtemp(prefix="espritz_res_"))
    try:
        for r in records:
            fasta = workdir / f"{_safe_stem(r['accession'])}.fasta"
            fasta.write_text(
                f">{r['accession']}\n{str(r['sequence']).upper()}\n", encoding="utf-8"
            )
        result = subprocess.run(
            [str(wrapper), str(workdir), model, str(sw)],
            capture_output=True, text=True, timeout=300,
            env=contract_env(),
        )
        if result.returncode != 0:
            return {}
        out: dict[str, dict[str, object]] = {}
        for r in records:
            path = workdir / f"{_safe_stem(r['accession'])}.espritz"
            if not path.exists():
                continue
            states, scores = _parse_espritz_residues(path)
            out[r["accession"]] = {"label": states, "residue": scores}
        return out
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _seg_masked_map(seg_bin: Path, fasta_path: Path) -> dict[str, str]:
    result = subprocess.run(
        [str(seg_bin), str(fasta_path), "12", "2.2", "2.5", "-x"],
        capture_output=True, text=True, timeout=120,
        env=contract_env(),
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
    return masked


def seg_residues(records: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    """Per-residue SEG low-complexity ``label`` mask (1 = masked)."""
    seg_bin = find_tool_entry("seg")
    fasta_text = "".join(
        f">{r['accession']}\n{str(r['sequence']).upper()}\n" for r in records
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="seg_res_"
    ) as fh:
        fh.write(fasta_text)
        fasta_path = Path(fh.name)
    try:
        masked = _seg_masked_map(seg_bin, fasta_path)
        out: dict[str, dict[str, object]] = {}
        for r in records:
            seq = masked.get(r["accession"], "")
            label = [1.0 if c == "x" or c.islower() else 0.0 for c in seq]
            out[r["accession"]] = {"label": label}
        return out
    finally:
        fasta_path.unlink(missing_ok=True)


def deepcoil_residues(records: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    """Per-residue DeepCoil ``coiled-coil`` (raw_cc) + ``sharpen`` (cc) arrays."""
    if not records:
        return {}
    runner = find_tool_entry_or_none("deepcoil")
    if runner is None:
        return {}
    out: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory(prefix="deepcoil_res_") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "input.fasta"
        out_dir = workdir / "out"
        out_dir.mkdir()
        lines: list[str] = []
        id_map: dict[str, str] = {}
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
                [str(runner), "-i", str(fasta), "-out_path", str(out_dir), "-n_cpu", "4"],
                check=True, capture_output=True, text=True,
                timeout=max(300, 30 * len(records)),
                env=contract_env(),
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return {}
        for safe, orig in id_map.items():
            fpath = out_dir / f"{safe}.out"
            if not fpath.exists():
                continue
            try:
                frame = pd.read_csv(fpath, sep="\t")
                out[orig] = {
                    "coiled-coil": frame["raw_cc"].astype(float).tolist(),
                    "sharpen": frame["cc"].astype(float).tolist(),
                }
            except Exception:
                continue
    return out


# ---------------------------------------------------------------------------
# Lookups (DeepPhase table, PhosphoSitePlus)
# ---------------------------------------------------------------------------
def lookup_deepphase(accessions: list[str]) -> dict[str, float]:
    """Look up DeepPhase scores from the bundled TSV or the supplement xlsx.

    DeepPhase has no inference tool — only a reference table. Prefers the
    package-bundled derived TSV; falls back to the table S3 xlsx. Accessions
    not in the table return nothing (the caller imputes with a warning).
    """
    table = deepphase_table_path()
    if not table.exists():
        warnings.warn(
            (
                f"DeepPhase table missing at {table}. "
                "Human-mode predictions (hSaPS / hPdPS) will use median "
                "imputation for the DeepPhase column. See the README "
                "section 'Optional human-feature data' for the download."
            ),
            PhaSePredMissingDataWarning,
            stacklevel=2,
        )
        return {}
    if table.suffix == ".tsv":
        dp = pd.read_csv(table, sep="\t")
        id_column = "Swiss-Prot_ID"
    else:
        dp = pd.read_excel(str(table), sheet_name="tableS3")
        id_column = "Swiss-Prot ID"
    dp_map: dict[str, float] = {}
    for _, row in dp.iterrows():
        swiss = str(row.get(id_column, "")).strip()
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


def compute_phos_freq(records: list[dict[str, str]]) -> dict[str, float]:
    """Compute phosphorylation-site frequency from PhosphoSitePlus.

    Returns per-accession site count / sequence length for the human set.
    """
    path = phosphosite_path()
    if not path.exists():
        warnings.warn(
            (
                f"PhosphoSitePlus dataset missing at {path}. "
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

    seq_rows: list[dict[str, object]] = [
        {"UniprotEntry": r["accession"], "sequence": r["sequence"]} for r in records
    ]
    features = compute_phosphosite_frequencies(seq_rows, path)
    hits = {f.accession: f.value for f in features}
    requested = {r["accession"] for r in records}
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