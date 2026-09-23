"""Command-line interface for the catgranule package.

``catgranule predict --fasta <in.fa> [--out json|csv] [--raw] [--route <r>]``
reads FASTA records and writes per-protein scores to stdout.

Weight-selection precedence: ``--weights <path>`` > ``$CATGRANULE_WEIGHTS``
> ``--route`` > packaged default (distilled).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Annotated

import typer

from catgranule import __version__ as _package_version
from catgranule.scoring import CatGranuleInputError, score_sequence
from catgranule.weights import (
    ROUTE_ARTIFACTS,
    CatGranuleWeights,
    CatGranuleWeightsError,
    load_weights,
)

app = typer.Typer(no_args_is_help=True)

RouteOption = Annotated[
    str | None,
    typer.Option(
        "--route",
        "-r",
        help="Weights artifact route: distilled (default), paper, legacy-unaudited.",
    ),
]
WeightsOption = Annotated[
    Path | None,
    typer.Option("--weights", help="Path to a weights artifact JSON (overrides env/route)."),
]


def _read_fasta(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    lines: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(lines)))
            header = line[1:].split(maxsplit=1)[0]
            lines = []
        elif header is not None and line:
            lines.append(line)
    if header is not None:
        records.append((header, "".join(lines)))
    return records


def _resolve_weights(route: str | None, weights_path: Path | None) -> CatGranuleWeights:
    if route is not None and route not in ROUTE_ARTIFACTS:
        raise typer.BadParameter(
            f"--route must be one of {sorted(ROUTE_ARTIFACTS)!r}, got {route!r}"
        )
    try:
        return load_weights(route=route, path=weights_path)
    except CatGranuleWeightsError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _output_json(
    records: list[tuple[str, str]], *, raw_mode: bool, weights: CatGranuleWeights
) -> None:
    rows: list[dict[str, object]] = []
    for accession, sequence in records:
        try:
            scored = score_sequence(sequence, normalized=not raw_mode, weights=weights)
            rows.append({"accession": accession, **scored})
        except CatGranuleInputError as exc:
            rows.append(
                {
                    "accession": accession,
                    "single": None,
                    "residue": None,
                    "granule_strength": None,
                    "error": str(exc),
                }
            )
    sys.stdout.write(json.dumps(rows, indent=2) + "\n")


def _output_csv(
    records: list[tuple[str, str]], *, raw_mode: bool, weights: CatGranuleWeights
) -> None:
    import csv as _csv

    writer = _csv.writer(sys.stdout)
    writer.writerow(["accession", "single", "granule_strength", "residue_length"])
    for accession, sequence in records:
        try:
            scored = score_sequence(sequence, normalized=not raw_mode, weights=weights)
            residue = scored["residue"]
            writer.writerow(
                [
                    accession,
                    scored["single"],
                    scored["granule_strength"],
                    len(residue) if isinstance(residue, (list, tuple)) else "",
                ]
            )
        except CatGranuleInputError:
            writer.writerow([accession, "", "", ""])


@app.command("predict")
def predict(
    fasta: Path = typer.Option(..., "--fasta", "-f", help="Input FASTA file."),
    out: str = typer.Option("json", "--out", help="Output format: json or csv."),
    raw: bool = typer.Option(
        False,
        "--raw",
        help="Emit raw Equation-2 values instead of Z-normalized (distilled "
        "route: the un-aligned formula profile; single unchanged).",
    ),
    route: RouteOption = None,
    weights: WeightsOption = None,
) -> None:
    """Score every protein in a FASTA file (single + residue profile)."""
    text = fasta.read_text(encoding="utf-8")
    records = _read_fasta(text)
    if not records:
        raise typer.BadParameter(f"No FASTA records found in {fasta}")
    resolved = _resolve_weights(route, weights)
    fmt = out.strip().lower()
    if fmt == "json":
        _output_json(records, raw_mode=raw, weights=resolved)
    elif fmt == "csv":
        _output_csv(records, raw_mode=raw, weights=resolved)
    else:
        raise typer.BadParameter("--out must be 'json' or 'csv'")


@app.command("check")
def check(
    route: RouteOption = None,
    weights: WeightsOption = None,
) -> None:
    """Show the loaded weights artifact and its provenance."""
    resolved = _resolve_weights(route, weights)
    provenance = resolved.provenance
    lines = [
        f"artifact_version: {resolved.version}  (the shipped weights artifact's own stamp)",
        f"package_version: {_package_version}",
        f"mode: {resolved.mode}",
        f"provenance.route: {provenance.route}",
        f"provenance.source: {provenance.source}",
        f"provenance.date: {provenance.date}",
        f"provenance.gate: {json.dumps(provenance.gate) if provenance.gate else None}",
        f"scales: {len(resolved.scales)} residues",
        f"normalization: Z({resolved.normalization_mean:.6f}, {resolved.normalization_std:.6f})",
    ]
    if resolved.model_path is not None:
        lines.append(f"model: {resolved.model_path}")
    typer.echo("\n".join(lines))


def main() -> None:
    app()