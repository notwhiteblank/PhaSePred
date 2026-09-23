from __future__ import annotations

import csv
import math
from io import StringIO
from pathlib import Path
from typing import Annotated

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import requests
import typer
from requests import HTTPError, RequestException
from sklearn.metrics import average_precision_score, roc_auc_score

from phasepred import __version__ as _phasepred_version
from phasepred import data
from phasepred.features import FeatureSchemaError, compute_native_features
from phasepred.legacy_features import (
    compute_lcr_features,
    compute_phosphosite_frequencies,
    compute_plaac_features,
    compute_pscore_features,
    load_deepphase_scores,
    native_feature_frame,
)
from phasepred.models import (
    load_model_artifact,
    predict_with_artifact,
    train_model_artifact,
)
from phasepred.paper import load_paper_features, load_paper_split
from phasepred.sequences import (
    UNIPROT_FASTA_BASE_URL,
    ProteinSequence,
    read_fasta_records,
    resolve_accessions,
    sequence_records_to_frame,
)

app = typer.Typer(no_args_is_help=True)


@app.callback(invoke_without_command=True)
def _root(
    ctx: typer.Context,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            help="Print the phasepred version and exit.",
            is_eager=True,
        ),
    ] = False,
) -> None:
    """PhaSePred — paper-faithful phase-separation predictor (Chen et al. 2022)."""
    # Single source of truth: this reports ``phasepred.__version__``, which
    # ``scripts/check_docs_consistency.py`` class (a) already pins to
    # ``pyproject.toml``/``pixi.toml``. It is deliberately not a fourth version
    # literal.
    if version:
        typer.echo(_phasepred_version)
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()

# Optional human-mode data (PhosphoSitePlus / DeepPhase) is resolved inside the
# command body through phasepred.tools, never as a module-level default. A
# frozen data_path(...) constant is evaluated at import time and points under
# lib/python3.*/ in a wheel install, so Click's exists=True would reject the
# command even when the data is correctly installed in the user data tiers.


TEST_AUC_ROWS = [
    {
        "dataset": "SaPS-test",
        "PScore": 0.80,
        "PLAAC": 0.82,
        "catGRANULE": 0.80,
        "FuzDrop": 0.81,
        "PhaSePred": 0.86,
        "PhaSePred remake": 0.8202521101294108,
    },
    {
        "dataset": "PdPS-test",
        "PScore": 0.61,
        "PLAAC": 0.62,
        "catGRANULE": 0.65,
        "FuzDrop": 0.59,
        "PhaSePred": 0.66,
        "PhaSePred remake": 0.7111709980362504,
    },
    {
        "dataset": "hSaPS-test",
        "PScore": 0.79,
        "PLAAC": 0.79,
        "catGRANULE": 0.77,
        "FuzDrop": 0.74,
        "PhaSePred": 0.85,
        "PhaSePred remake": 0.8694098883572569,
    },
    {
        "dataset": "hPdPS-test",
        "PScore": 0.58,
        "PLAAC": 0.58,
        "catGRANULE": 0.66,
        "FuzDrop": 0.55,
        "PhaSePred": 0.78,
        "PhaSePred remake": 0.8379463307776561,
    },
]


@app.command("prepare-split")
def prepare_split(
    s2: Annotated[Path, typer.Option("--s2", exists=True, readable=True)],
    s3: Annotated[Path, typer.Option("--s3", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    split = load_paper_split(s2, s3)
    output.parent.mkdir(parents=True, exist_ok=True)
    split.all_rows.to_csv(output, index=False)
    typer.echo(f"Wrote {len(split.all_rows)} split rows to {output}")


@app.command("features-from-fasta")
def features_from_fasta(
    input_path: Annotated[Path, typer.Option("--input", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    """Compute v2022 sequence features from a FASTA file.

    Emits native features (Hydropathy, FCR) plus the ESpritz-domain IDR
    fraction and the SEG-domain LCR fraction. Missing external tools degrade
    the corresponding column to NaN with a warning (v2022 definitions).
    """
    records = read_fasta_records(input_path)
    rows = []
    idr_rows = _compute_espritz_idr_fraction(records)
    lcr_rows = _compute_lcr_fraction(records)
    for record in records:
        row = {
            "UniprotEntry": record.accession,
            "sequence": record.sequence,
            **compute_native_features(record.sequence),
        }
        row["espritz-idr-fraction"] = idr_rows.get(record.accession, float("nan"))
        row["lcr-fraction"] = lcr_rows.get(record.accession, float("nan"))
        rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    typer.echo(f"Wrote {len(rows)} v2022 feature rows to {output}")


def _compute_espritz_idr_fraction(records) -> dict[str, float]:
    """IDR fraction (ESpritz DisProt, 5% FPR) with NaN degradation.

    Emits a warning on stderr so the NaN degradation is always visible.
    """
    try:
        from phasepred.tools import run_espritz

        return run_espritz([{"accession": r.accession, "sequence": r.sequence} for r in records])
    except Exception as exc:
        typer.echo(
            f"WARNING: ESpritz unavailable; 'espritz-idr-fraction' set to NaN "
            f"({exc}).",
            err=True,
        )
        return {}


def _compute_lcr_fraction(records) -> dict[str, float]:
    """LCR fraction (SEG default parameters) with NaN degradation."""
    try:
        from phasepred.tools import run_seg

        return run_seg([{"accession": r.accession, "sequence": r.sequence} for r in records])
    except Exception as exc:
        typer.echo(
            f"WARNING: SEG unavailable; 'lcr-fraction' set to NaN ({exc}).",
            err=True,
        )
        return {}


@app.command("sequences-from-ids")
def sequences_from_ids(
    input_path: Annotated[Path, typer.Option("--input", exists=True, readable=True)],
    cache: Annotated[Path, typer.Option("--cache")],
    output: Annotated[Path, typer.Option("--output")],
    allow_network: Annotated[bool, typer.Option("--allow-network/--no-network")] = True,
) -> None:
    frame = _read_table(input_path)
    if "UniprotEntry" not in frame.columns:
        raise typer.BadParameter("Input table must contain a UniprotEntry column")
    accessions = list(dict.fromkeys(frame["UniprotEntry"].dropna().astype(str).str.strip()))
    records = resolve_accessions(
        accessions,
        cache_path=cache,
        allow_network=allow_network,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    sequence_records_to_frame(records).to_csv(output, index=False)
    typer.echo(f"Wrote {len(records)} sequence records to {output}")


@app.command("features-from-paper")
def features_from_paper(
    s2: Annotated[Path, typer.Option("--s2", exists=True, readable=True)],
    s3: Annotated[Path, typer.Option("--s3", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    features = load_paper_features(s2, s3)
    output.parent.mkdir(parents=True, exist_ok=True)
    features.all_rows.to_csv(output, index=False)
    typer.echo(f"Wrote {len(features.all_rows)} paper feature rows to {output}")


@app.command("features-from-split")
def features_from_split(
    split_table: Annotated[Path, typer.Option("--split-table", exists=True, readable=True)],
    sequences: Annotated[Path, typer.Option("--sequences", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    split = _read_table(split_table)
    sequence_frame = _read_table(sequences)
    merged = split.merge(sequence_frame, on="UniprotEntry", how="left")
    if merged["sequence"].isna().any():
        missing = merged.loc[merged["sequence"].isna(), "UniprotEntry"].tolist()
        raise typer.BadParameter(f"Missing sequences for accessions: {', '.join(missing)}")
    rows = []
    for row in merged.itertuples(index=False):
        features = {
            "task": row.task,
            "split": row.split,
            "label": row.label,
            "UniprotEntry": row.UniprotEntry,
            "sequence": row.sequence,
            **compute_native_features(row.sequence),
        }
        rows.append(features)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    typer.echo(f"Wrote {len(rows)} split feature rows to {output}")


def _resolve_human_data_paths(
    phosphosite: Path | None, deepphase: Path | None
) -> tuple[Path, Path]:
    """Resolve the two human-mode data paths at call time (never at import).

    Both options default to ``None``; an explicit flag wins. When omitted,
    ``tools.phosphosite_path`` / ``tools.deepphase_table_path`` re-read the
    environment and walk the three-tier cascade (PHASEPRED_DATA_ROOT >
    XDG/home user data > checkout), so a wheel install with tier-1/tier-2 data
    works. A missing file raises a message naming the cascade and the installer
    instead of Click's generic "does not exist".
    """
    from phasepred import tools

    resolved_phosphosite = phosphosite or tools.phosphosite_path()
    resolved_deepphase = deepphase or tools.deepphase_table_path()
    if not resolved_phosphosite.is_file():
        raise typer.BadParameter(
            "PhosphoSitePlus dataset not found. Resolution order: "
            "$PHASEPRED_DATA_ROOT/phosphositeplus/, "
            "$XDG_DATA_HOME/phasepred/phosphositeplus/, then the checkout's "
            "data/raw/external/phosphositeplus/. Install it with: "
            "bash tools/PhosphoSitePlus/install.sh"
        )
    if not resolved_deepphase.is_file():
        raise typer.BadParameter(
            "DeepPhase score table not found. Resolution order: the wheel-bundled "
            "data/deepphase/deepphase_scores.tsv, then "
            "$PHASEPRED_DATA_ROOT/deepphase/extracted/tableS3.xlsx, "
            "$XDG_DATA_HOME/phasepred/deepphase/..., then the checkout's "
            "data/raw/external/deepphase/extracted/. Reinstall the phasepred "
            "wheel or re-derive the table with scripts/derive_deepphase_tsv.py"
        )
    return resolved_phosphosite, resolved_deepphase


@app.command("features-recomputed")
def features_recomputed(
    sequences: Annotated[Path, typer.Option("--sequences", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
    espritz_cache: Annotated[Path | None, typer.Option("--espritz-cache")] = None,
    include_human: Annotated[bool, typer.Option("--include-human/--base-only")] = False,
    phosphosite: Annotated[
        Path | None,
        typer.Option(
            "--phosphosite",
            help="PhosphoSitePlus dataset; defaults to the user-data cascade.",
        ),
    ] = None,
    deepphase: Annotated[
        Path | None,
        typer.Option(
            "--deepphase",
            help="DeepPhase score table; defaults to the wheel-bundled TSV.",
        ),
    ] = None,
    catgranule_features: Annotated[
        Path | None,
        typer.Option("--catgranule-features", exists=True, readable=True),
    ] = None,
    catgranule_source: Annotated[
        str,
        typer.Option(
            "--catgranule-source",
            help="catGRANULE source: 'package' (catgranule pkg) or 'csv'.",
        ),
    ] = "package",
    deepcoil_features: Annotated[
        Path | None,
        typer.Option("--deepcoil-features", exists=True, readable=True),
    ] = None,
    pscore_features: Annotated[
        Path | None,
        typer.Option("--pscore-features", exists=True, readable=True),
    ] = None,
    plaac_features: Annotated[
        Path | None,
        typer.Option("--plaac-features", exists=True, readable=True),
    ] = None,
) -> None:
    human_phosphosite: Path | None = phosphosite
    human_deepphase: Path | None = deepphase
    if include_human:
        human_phosphosite, human_deepphase = _resolve_human_data_paths(phosphosite, deepphase)
    sequence_frame = _read_table(sequences)
    if "UniprotEntry" not in sequence_frame.columns or "sequence" not in sequence_frame.columns:
        raise typer.BadParameter("Sequence table must contain UniprotEntry and sequence columns")
    rows = sequence_frame[["UniprotEntry", "sequence"]].drop_duplicates("UniprotEntry")
    sequence_rows = rows.to_dict(orient="records")
    feature_frame = native_feature_frame(sequence_rows)

    lcr = pd.DataFrame([f.to_json() for f in compute_lcr_features(sequence_rows)])
    if pscore_features is not None:
        pscore = _read_table(pscore_features)
        if "UniprotEntry" not in pscore.columns or "PScore" not in pscore.columns:
            raise typer.BadParameter(
                "--pscore-features must contain UniprotEntry and PScore columns"
            )
    else:
        pscore = pd.DataFrame([f.to_json() for f in compute_pscore_features(sequence_rows)])
    if plaac_features is not None:
        plaac = _read_table(plaac_features)
        if "UniprotEntry" not in plaac.columns or "PLAAC" not in plaac.columns:
            raise typer.BadParameter("--plaac-features must contain UniprotEntry and PLAAC columns")
    else:
        plaac = pd.DataFrame([f.to_json() for f in compute_plaac_features(sequence_rows)])
    feature_frame = _merge_scalar_feature(feature_frame, lcr, "LCR")
    feature_frame = _merge_scalar_feature(feature_frame, pscore, "PScore")
    feature_frame = _merge_scalar_feature(feature_frame, plaac, "PLAAC")

    if espritz_cache is None:
        raise typer.BadParameter("--espritz-cache is required (ESpritz is the only IDR source)")
    from phasepred.espritz import ESpritzFeatureCache

    idr = pd.DataFrame(
        [f.to_json() for f in ESpritzFeatureCache(espritz_cache)._read_all().values()]
    ).rename(columns={"accession": "UniprotEntry", "espritz_idr": "IDR"})
    feature_frame = feature_frame.merge(idr[["UniprotEntry", "IDR"]], on="UniprotEntry", how="left")

    catgranule_source_normalized = catgranule_source.strip().lower()
    if catgranule_source_normalized == "csv":
        if catgranule_features is None:
            raise typer.BadParameter(
                "--catgranule-features is required when --catgranule-source csv"
            )
        catgranule_frame = _read_table(catgranule_features)
    elif catgranule_source_normalized == "package":
        from catgranule import CatGranuleInputError, score_sequence

        catgranule_frame = _compute_catgranule_package_features(
            sequence_rows, score_sequence, CatGranuleInputError
        )
    else:
        raise typer.BadParameter("--catgranule-source must be 'package' or 'csv'")
    if deepcoil_features is None:
        raise typer.BadParameter(
            "--deepcoil-features is required for full 8/10 feature recomputation; "
            "DeepCoil is not installed in the project Python 3.12 environment."
        )
    deepcoil_frame = _read_table(deepcoil_features)
    for label, frame, column in (
        ("catgranule", catgranule_frame, "catGRANULE"),
        ("deepcoil", deepcoil_frame, "DeepCoil"),
    ):
        missing = {"UniprotEntry", column} - set(frame.columns)
        if missing:
            raise typer.BadParameter(
                f"{label} feature table missing required columns: {', '.join(sorted(missing))}"
            )
    feature_frame = feature_frame.merge(
        catgranule_frame[["UniprotEntry", "catGRANULE"]],
        on="UniprotEntry",
        how="left",
    )
    feature_frame = feature_frame.merge(
        deepcoil_frame[["UniprotEntry", "DeepCoil"]],
        on="UniprotEntry",
        how="left",
    )

    if include_human:
        assert human_phosphosite is not None and human_deepphase is not None
        phospho = pd.DataFrame(
            [
                f.to_json()
                for f in compute_phosphosite_frequencies(sequence_rows, human_phosphosite)
            ]
        )
        feature_frame = _merge_scalar_feature(feature_frame, phospho, "Phos freq")
        feature_frame = feature_frame.merge(
            load_deepphase_scores(human_deepphase),
            on="UniprotEntry",
            how="left",
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    feature_frame.to_csv(output, index=False)
    typer.echo(f"Wrote {len(feature_frame)} recomputed feature rows to {output}")


@app.command("train")
def train(
    split_table: Annotated[Path, typer.Option("--split-table", exists=True, readable=True)],
    features: Annotated[Path, typer.Option("--features", exists=True, readable=True)],
    task: Annotated[str, typer.Option("--task")],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    split = _read_table(split_table)
    feature_frame = _read_table(features)
    task_split = split.query("task == @task and split == 'train'")
    if task_split.empty:
        raise typer.BadParameter(f"No train rows found for task {task}")
    merge_keys = ["UniprotEntry"]
    if "task" in feature_frame.columns:
        merge_keys.append("task")
    train_frame = task_split.merge(feature_frame, on=merge_keys, how="left")
    artifact = train_model_artifact(train_frame, label_column="label", task=task)
    artifact.dump(output)
    typer.echo(f"Wrote {task} model to {output}")


@app.command("predict-features")
def predict_features(
    model: Annotated[Path, typer.Option("--model", exists=True, readable=True)],
    features: Annotated[Path, typer.Option("--features", exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    artifact = load_model_artifact(model)
    frame = _read_table(features)
    predictions = predict_with_artifact(artifact, frame)
    output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output, index=False)
    typer.echo(f"Wrote {len(predictions)} predictions to {output}")


@app.command("predict")
def predict(
    mode: Annotated[str, typer.Option("--mode", help="SaPS, PdPS, hSaPS, or hPdPS")],
    output: Annotated[Path, typer.Option("--output")],
    fasta: Annotated[
        Path | None,
        typer.Option("--fasta", help="FASTA file with protein sequences."),
    ] = None,
    ids: Annotated[
        str | None,
        typer.Option(
            "--ids",
            help="UniProt accessions. Three accepted forms: a comma-separated "
                 "inline list ('P35637,Q9Y2W1'), a path to a file with one ID "
                 "per line, or '-' to read IDs from stdin. Sequences are fetched "
                 "from rest.uniprot.org and cached locally.",
        ),
    ] = None,
    cache: Annotated[
        Path | None,
        typer.Option(
            "--cache",
            help="JSONL cache for UniProt fetches (default: data/interim/uniprot_cache.jsonl).",
        ),
    ] = None,
    models_dir: Annotated[
        Path | None,
        typer.Option(
            "--models-dir",
            help="Override the default Product A model dir entirely. Expects "
                 "<dir>/<mode>/8f_model_*.joblib.",
        ),
    ] = None,
) -> None:
    """Predict phase separation scores for proteins.

    Provide proteins via `--fasta`, `--ids`, or both. `--ids` accepts:
      * a comma-separated inline list: `--ids P35637,Q9Y2W1`
      * a file with one accession per line: `--ids my_ids.txt`
      * standard input: `cat my_ids.txt | phasepred predict --ids - ...`

    Examples:
        phasepred predict --fasta my_proteins.fasta --mode SaPS --output scores.csv
        phasepred predict --ids "P35637,Q9Y2W1" --mode hSaPS --output scores.csv
    """
    from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS
    from phasepred.models import model_dir_is_stale
    from phasepred.predictor import (
        DEFAULT_MODELS_DIR,
        compute_all_features,
        load_ensemble,
        load_single_model,
        predict_scores,
    )
    from phasepred.sequences import resolve_accessions
    from phasepred.tool_paths import (
        PhaSePredToolNotFound,
        assert_predict_requirements,
    )

    if fasta is None and ids is None:
        raise typer.BadParameter("Provide at least one of --fasta or --ids.")

    mode = mode.strip()
    if mode not in ("SaPS", "PdPS", "hSaPS", "hPdPS"):
        raise typer.BadParameter(
            f"Unknown mode: {mode}. Must be SaPS, PdPS, hSaPS, or hPdPS."
        )

    chosen_models_dir = models_dir if models_dir is not None else DEFAULT_MODELS_DIR

    # Preflight: every external tool that `predict` needs must resolve up-front.
    # If anything is MISSING the user gets a clear error before we read any input.
    try:
        assert_predict_requirements()
    except PhaSePredToolNotFound as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1) from e

    records: list[dict[str, str]] = []
    seen: set[str] = set()

    # --- FASTA input ---
    if fasta is not None:
        if not fasta.exists() or not fasta.is_file():
            raise typer.BadParameter(f"--fasta path not found or not a file: {fasta}")
        fasta_text = Path(fasta).read_text(encoding="utf-8")
        fasta_records = list(_parse_fasta_records(fasta_text))
        if not fasta_records:
            raise typer.BadParameter(f"No FASTA records found in {fasta}")
        for acc, seq, _hdr in fasta_records:
            if acc in seen:
                continue
            records.append({"accession": acc, "sequence": seq})
            seen.add(acc)
        typer.echo(f"Read {len(fasta_records)} record(s) from {fasta}")

    # --- IDs input ---
    if ids is not None:
        accessions = _resolve_ids_argument(ids)
        accessions = [a for a in accessions if a and a not in seen]
        if accessions:
            cache_path = cache if cache is not None else data.data_path(
                "data", "interim", "uniprot_cache.jsonl"
            )
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            typer.echo(f"Resolving {len(accessions)} UniProt ID(s) via {cache_path}...")
            try:
                fetched = resolve_accessions(accessions, cache_path=cache_path, allow_network=True)
            except Exception as e:
                typer.echo(f"UniProt resolution failed: {e}", err=True)
                raise typer.Exit(code=1) from e
            for rec in fetched:
                if rec.accession in seen:
                    continue
                records.append({"accession": rec.accession, "sequence": rec.sequence})
                seen.add(rec.accession)
            typer.echo(f"Resolved {len(fetched)} ID(s) to sequences")

    if not records:
        raise typer.BadParameter(
            "No protein records to score after merging --fasta and --ids inputs."
        )

    # Compute features
    is_human = mode.startswith("h")
    typer.echo(f"Computing features for {len(records)} protein(s), mode={mode}...")
    feature_df = compute_all_features(records, include_human=is_human)

    # Load model(s)
    model_path_or_dir = chosen_models_dir / f"{mode}.joblib"
    if model_path_or_dir.exists():
        typer.echo(f"Loading model from {model_path_or_dir}")
        model = load_single_model(model_path_or_dir)
        models = [model]
        stale = True
        stale_reason = (
            f"model artifact {model_path_or_dir} carries no feature-definition "
            "version tag (legacy or externally trained model)"
        )
    else:
        model_dir = chosen_models_dir / mode
        if not model_dir.exists():
            raise typer.BadParameter(
                f"Model not found at {model_path_or_dir} or {model_dir}.\n"
                f"Ensure the Product A models exist under {DEFAULT_MODELS_DIR}, "
                "or use --models-dir to override."
            )
        typer.echo(f"Loading ensemble from {model_dir}")
        models = load_ensemble(model_dir)
        stale, stale_reason = model_dir_is_stale(model_dir)

    if stale:
        typer.echo(
            "\n".join(
                [
                    "WARNING: model features are STALE.",
                    f"  {stale_reason}.",
                    "  Scores are produced with feature definitions that differ "
                    "from the current code (feature_definitions_version v2022).",
                    "  The output CSV carries a 'stale' column so downstream "
                    "consumers can detect this.",
                ]
            ),
            err=True,
        )

    # Predict
    columns = HUMAN_FEATURE_COLUMNS if is_human else BASE_FEATURE_COLUMNS
    predictions = predict_scores(models, feature_df, list(columns))

    # Merge features for informative output
    output_df = predictions.merge(feature_df, on="UniprotEntry", how="left")
    output_df["stale"] = stale
    output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output, index=False)

    for _, row in output_df.iterrows():
        typer.echo(f"  {row['UniprotEntry']}: score={row['score']:.4f}")
    typer.echo(f"Wrote {len(output_df)} predictions to {output}")


def _resolve_ids_argument(raw: str) -> list[str]:
    """Resolve the --ids argument to a list of accessions.

    Accepts:
      * '-': read whitespace-separated tokens from stdin
      * existing file path: one accession per line (blank lines / # comments allowed)
      * otherwise: comma-separated inline string
    Whitespace is stripped; duplicates are deduplicated preserving order.
    """
    import sys as _sys

    raw = raw.strip()
    if raw == "-":
        tokens = _sys.stdin.read().split()
    else:
        candidate = Path(raw)
        if candidate.is_file():
            tokens = []
            for line in candidate.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                tokens.append(line)
        else:
            tokens = [tok.strip() for tok in raw.split(",")]
    seen: dict[str, None] = {}
    for tok in tokens:
        if tok and tok not in seen:
            seen[tok] = None
    return list(seen.keys())


@app.command("check-tools")
def check_tools(
    strict: Annotated[
        bool,
        typer.Option(
            "--strict/--no-strict",
            help="Exit with code 1 if any required-for-predict tool is MISSING.",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Emit a JSON array on stdout; all diagnostics go to stderr.",
        ),
    ] = False,
) -> None:
    """Report the install status of every external feature tool.

    Each tool is a contract package under tools/<TOOL>/ (run +
    manifest.toml + install.sh + README.md). Resolution cascade per tool:
    PHASEPRED_<TOOL>_DIR, the user data directory, the vendored tools/<TOOL>/
    package, PATH (manifest path_executables), then missing. Usability is
    verified with each package's ./run --check. Run install.sh after
    installing. ``--json`` keeps stdout to a single parseable JSON array
    (Ruling 4); anything a probe prints goes to stderr.
    """
    import contextlib
    import sys as _sys

    from phasepred.tool_paths import (
        check_all_tools,
        format_status_table,
        statuses_to_json,
    )

    if json_output:
        # Ruling 4: only the JSON may reach stdout. Redirect any probe's
        # stdout (diagnostics/progress) to stderr while statuses are computed.
        with contextlib.redirect_stdout(_sys.stderr):
            statuses = check_all_tools()
        typer.echo(statuses_to_json(statuses))
    else:
        statuses = check_all_tools()
        typer.echo(format_status_table(statuses))

    if strict:
        missing_required = [s for s in statuses if not s.ok and s.required_for_predict]
        if missing_required:
            raise typer.Exit(code=1)


@app.command("tool-features")
def tool_features(
    fasta: Annotated[
        Path | None,
        typer.Option("--input", help="FASTA file with proteins to decompose."),
    ] = None,
    ids: Annotated[
        str | None,
        typer.Option(
            "--ids",
            help="UniProt accessions: inline comma list, a file (one per line), "
            "or '-' for stdin.",
        ),
    ] = None,
    tools_arg: Annotated[
        str,
        typer.Option(
            "--tools",
            help="Comma-separated tool list, or 'all'. "
            "Vocabulary: espritz,seg,pscore,plaac,catgranule,deepcoil,"
            "deepphase,phos,hydropathy,charge.",
        ),
    ] = "all",
    output: Annotated[Path, typer.Option("--output", help="CSV summary output.")] = Path(
        "tool-features.csv"
    ),
    residue_json: Annotated[
        Path | None,
        typer.Option(
            "--residue-json",
            help="Directory for per-protein per-tool residue-level JSON.",
        ),
    ] = None,
    cache: Annotated[Path | None, typer.Option("--cache")] = None,
    allow_network: Annotated[bool, typer.Option("--allow-network/--no-network")] = True,
) -> None:
    """Decompose each protein's PS profile tool-by-tool.

    Summary values per tool mirror the web ``tools.<x>`` payloads; a single
    failing tool degrades to NaN (never aborts). ``--residue-json`` adds
    per-residue arrays for espritz/seg/catgranule/deepcoil/charge.
    """
    from phasepred.sequences import resolve_accessions
    from phasepred.tool_features import (
        TOOL_NAMES,
        compute_tool_features,
        write_residue_json,
    )

    if fasta is None and ids is None:
        raise typer.BadParameter("Provide at least one of --input or --ids.")

    requested = (
        list(TOOL_NAMES)
        if tools_arg.strip().lower() == "all"
        else [t.strip() for t in tools_arg.split(",") if t.strip()]
    )
    unknown = [t for t in requested if t not in TOOL_NAMES]
    if unknown:
        raise typer.BadParameter(
            f"Unknown tool(s): {', '.join(unknown)}. Valid: {', '.join(TOOL_NAMES)}"
        )

    records: list[dict[str, str]] = []
    seen: set[str] = set()

    if fasta is not None:
        if not fasta.exists():
            raise typer.BadParameter(f"--input path not found: {fasta}")
        for acc, seq, _hdr in _parse_fasta_records(fasta.read_text(encoding="utf-8")):
            if acc in seen:
                continue
            records.append({"accession": acc, "sequence": seq})
            seen.add(acc)

    if ids is not None:
        accessions = [a for a in _resolve_ids_argument(ids) if a and a not in seen]
        if accessions:
            cache_path = cache if cache is not None else data.data_path(
                "data", "interim", "uniprot_cache.jsonl"
            )
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            fetched = resolve_accessions(
                accessions, cache_path=cache_path, allow_network=allow_network
            )
            for rec in fetched:
                if rec.accession in seen:
                    continue
                records.append({"accession": rec.accession, "sequence": rec.sequence})
                seen.add(rec.accession)

    if not records:
        raise typer.BadParameter("No proteins to decompose after merging inputs.")

    typer.echo(f"Decomposing {len(records)} protein(s) across {len(requested)} tool(s)...")
    summary, residue_map = compute_tool_features(
        records, requested, residue=residue_json is not None
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, index=False)
    typer.echo(f"Wrote {len(summary)} summary row(s) to {output}")

    if residue_json is not None:
        write_residue_json(residue_map, residue_json)
        typer.echo(f"Wrote residue JSON to {residue_json}/")


def _read_table(path: Path) -> pd.DataFrame:
    text = Path(path).read_text(encoding="utf-8")
    sample = next((line for line in text.splitlines() if line.strip()), "")
    if not sample:
        return pd.DataFrame()
    delimiter = ","
    if "\t" in sample and "," not in sample:
        delimiter = "\t"
    else:
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t;|,")
            delimiter = dialect.delimiter
        except csv.Error:
            pass
    return pd.read_csv(StringIO(text), sep=delimiter)


def _frame_from_rows(
    rows: list[dict[str, float | str]],
    columns: list[str],
    *,
    index_column: str,
) -> pd.DataFrame:
    frame = pd.DataFrame(rows).set_index(index_column)
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise typer.BadParameter(f"Missing AUC columns: {', '.join(missing)}")
    return frame[columns]


def _compute_catgranule_package_features(
    sequence_rows: list[dict[str, object]],
    score_sequence_fn,
    error_cls,
) -> pd.DataFrame:
    rows = []
    for record in sequence_rows:
        accession = str(record["UniprotEntry"])
        sequence = str(record["sequence"])
        try:
            score = float(score_sequence_fn(sequence)["single"])
            rows.append({"UniprotEntry": accession, "catGRANULE": score})
        except error_cls:
            rows.append({"UniprotEntry": accession})
    return pd.DataFrame(rows)


def _merge_scalar_feature(
    frame: pd.DataFrame,
    scalar_frame: pd.DataFrame,
    column: str,
) -> pd.DataFrame:
    if scalar_frame.empty:
        return frame.assign(**{column: np.nan})
    values = scalar_frame.rename(columns={"accession": "UniprotEntry", "value": column})
    return frame.merge(values[["UniprotEntry", column]], on="UniprotEntry", how="left")


def _resolve_accessions_for_comparison(
    accessions: list[str],
    *,
    cache_path: Path,
    allow_network: bool,
    skip_missing: bool,
    missing_output: Path,
):
    if not skip_missing:
        return resolve_accessions(
            accessions,
            cache_path=cache_path,
            allow_network=allow_network,
        )
    from phasepred.sequences import SequenceCache

    cache = SequenceCache(cache_path)
    cached_records = cache._read_all()
    records: list[ProteinSequence] = []
    missing_rows: list[dict[str, str]] = []
    missing_accessions: list[str] = []
    seen: set[str] = set()
    for accession in accessions:
        normalized = accession.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        if normalized in cached_records:
            records.append(cached_records[normalized])
        else:
            missing_accessions.append(normalized)

    if not allow_network:
        missing_rows.extend(
            {"UniprotEntry": accession, "error": "not cached"} for accession in missing_accessions
        )
    elif missing_accessions:
        fetched_records, missing_rows = _fetch_uniprot_records_in_batches(
            missing_accessions,
            cache=cache,
            cached_records=cached_records,
        )
        for record in fetched_records:
            cached_records[record.accession] = record
            records.append(record)
    if missing_rows:
        missing_output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(missing_rows).to_csv(missing_output, index=False)
    return records


def _fetch_uniprot_records_in_batches(
    accessions: list[str],
    *,
    cache,
    cached_records: dict[str, ProteinSequence],
    batch_size: int = 50,
):
    from datetime import UTC, datetime

    fetched: list[ProteinSequence] = []
    missing_rows: list[dict[str, str]] = []
    session = requests.Session()
    for start in range(0, len(accessions), batch_size):
        batch = accessions[start : start + batch_size]
        try:
            response = session.get(
                f"{UNIPROT_FASTA_BASE_URL}/stream",
                params={
                    "compressed": "false",
                    "format": "fasta",
                    "query": " OR ".join(f"accession:{accession}" for accession in batch),
                },
                timeout=120,
            )
            response.raise_for_status()
        except HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else "unknown"
            missing_rows.extend(
                {"UniprotEntry": accession, "error": f"batch HTTP {status}"}
                for accession in batch
            )
            continue
        except RequestException as exc:
            missing_rows.extend(
                {"UniprotEntry": accession, "error": str(exc)} for accession in batch
            )
            continue

        returned = {}
        for accession, sequence, header in _parse_fasta_records(response.text):
            if accession in batch:
                fetched_record = ProteinSequence(
                    accession=accession,
                    sequence=sequence,
                    source="uniprot",
                    header=header,
                    retrieved_at=datetime.now(UTC).isoformat(),
                    endpoint=f"{UNIPROT_FASTA_BASE_URL}/stream",
                )
                returned[accession] = fetched_record
        for accession in batch:
            record: ProteinSequence | None = returned.get(accession)
            if record is None:
                missing_rows.append({"UniprotEntry": accession, "error": "not returned"})
                continue
            fetched.append(record)
            cached_records[record.accession] = record
        cache.write_all(cached_records)
    return fetched, missing_rows


def _parse_fasta_records(text: str):
    header: str | None = None
    sequence_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield _accession_from_fasta_header(header), "".join(sequence_lines), header
            header = line
            sequence_lines = []
        elif header is not None:
            sequence_lines.append(line.strip())
    if header is not None:
        yield _accession_from_fasta_header(header), "".join(sequence_lines), header


def _accession_from_fasta_header(header: str) -> str:
    identifier = header[1:].split(maxsplit=1)[0]
    parts = identifier.split("|")
    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return identifier


def _train_and_score_paper_features(
    features: pd.DataFrame,
    *,
    variant: str,
    model_dir: Path,
    prediction_dir: Path,
) -> list[dict[str, object]]:
    model_dir.mkdir(parents=True, exist_ok=True)
    prediction_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for task in ("SaPS", "PdPS", "hSaPS", "hPdPS"):
        train_frame = features.query("task == @task and split == 'train'").copy()
        test_frame = features.query("task == @task and split == 'test'").copy()
        if train_frame.empty:
            raise typer.BadParameter(f"No training rows found for task {task}")
        if test_frame.empty:
            raise typer.BadParameter(f"No test rows found for task {task}")
        artifact = train_model_artifact(train_frame, label_column="label", task=task)
        model_path = model_dir / f"{task}.joblib"
        artifact.dump(model_path)
        predictions = predict_with_artifact(artifact, test_frame)
        predictions.to_csv(prediction_dir / f"{task}_test_scores.csv", index=False)
        y_true = test_frame["label"].astype(int)
        y_score = predictions["score"]
        rows.append(
            {
                "variant": variant,
                "task": task,
                "model_path": str(model_path),
                "train_rows": len(train_frame),
                "test_rows": len(test_frame),
                "roc_auc": roc_auc_score(y_true, y_score)
                if y_true.nunique() > 1
                else math.nan,
                "average_precision": average_precision_score(y_true, y_score)
                if y_true.nunique() > 1
                else math.nan,
            }
        )
    return rows


def _native_features_from_sequences(sequence_frame: pd.DataFrame):
    rows = []
    errors = []
    for record in sequence_frame.itertuples(index=False):
        accession = str(record.UniprotEntry)
        row = {
            "UniprotEntry": accession,
            "sequence_hash": getattr(record, "sequence_hash", None),
        }
        try:
            row.update(compute_native_features(str(record.sequence)))
        except FeatureSchemaError as exc:
            errors.append({"UniprotEntry": accession, "error": str(exc)})
        rows.append(row)
    return rows, errors


def _train_and_score_feature_table(
    features: pd.DataFrame,
    *,
    variant: str,
    model_dir: Path,
    prediction_dir: Path,
    feature_columns: list[str],
) -> list[dict[str, object]]:
    model_dir.mkdir(parents=True, exist_ok=True)
    prediction_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for task in ("SaPS", "PdPS", "hSaPS", "hPdPS"):
        train_frame = features.query("task == @task and split == 'train'").copy()
        test_frame = features.query("task == @task and split == 'test'").copy()
        if train_frame.empty:
            raise typer.BadParameter(f"No training rows found for task {task}")
        if test_frame.empty:
            raise typer.BadParameter(f"No test rows found for task {task}")
        artifact = train_model_artifact(
            train_frame,
            label_column="label",
            task=task,
            feature_columns=feature_columns,
        )
        model_path = model_dir / f"{task}.joblib"
        artifact.dump(model_path)
        predictions = predict_with_artifact(artifact, test_frame)
        predictions.to_csv(prediction_dir / f"{task}_test_scores.csv", index=False)
        y_true = test_frame["label"].astype(int)
        y_score = predictions["score"]
        rows.append(
            {
                "variant": variant,
                "task": task,
                "model_path": str(model_path),
                "feature_columns": ";".join(feature_columns),
                "train_rows": len(train_frame),
                "test_rows": len(test_frame),
                "train_rows_with_missing_features": int(
                    train_frame[feature_columns].isna().any(axis=1).sum()
                ),
                "test_rows_with_missing_features": int(
                    test_frame[feature_columns].isna().any(axis=1).sum()
                ),
                "roc_auc": roc_auc_score(y_true, y_score)
                if y_true.nunique() > 1
                else math.nan,
                "average_precision": average_precision_score(y_true, y_score)
                if y_true.nunique() > 1
                else math.nan,
            }
        )
    return rows


def main() -> None:
    app()
