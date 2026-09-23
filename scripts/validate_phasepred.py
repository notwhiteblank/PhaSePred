#!/usr/bin/env python
"""One-command three-layer validation of the shipped PhaSePred artifacts.

Runs the AUC, leakage and artifact-consistency layers from
``phasepred.validation`` over the committed training tables and the published
models, writes a clock-free machine-readable report and returns a tiered exit
code (Ruling 1):

    0   no hard gate violated; non-blocking rows are reported as WARN
    1   a hard gate failed (accession-level leakage or artifact consistency);
        ``--strict-paper`` additionally promotes a paper-layer row beyond 0.01
        and ``--strict-leakage`` promotes the sequence-identity row
    2   input missing or usage error

A missing training table exits 2 with one actionable line rather than a
traceback (phases E6 and E7 branch on the code).

Usage:
    python scripts/validate_phasepred.py
    python scripts/validate_phasepred.py --sequences data/interim/uniprot_cache.jsonl
    python scripts/validate_phasepred.py --recompute-auc
    python scripts/validate_phasepred.py --strict-paper --json
    python scripts/validate_phasepred.py --sequences data/interim/uniprot_cache.jsonl \
        --strict-leakage
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from phasepred.training_data import TrainingDataSchemaError  # noqa: E402
from phasepred.validation import (  # noqa: E402
    ValidationError,
    exit_code_for,
    format_report_table,
    results_to_report,
    run_all,
)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metrics",
        type=Path,
        default=None,
        metavar="PATH",
        help="metrics.json to compare against the paper targets "
        "(default: products/A_paper_split_recomputed/metrics.json)",
    )
    parser.add_argument(
        "--models-dir",
        type=Path,
        default=None,
        metavar="PATH",
        help="models root holding the four task directories (default: bundled models root)",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        metavar="PATH",
        help="directory holding the training tables (default: the repository's data root)",
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        default=None,
        metavar="PATH",
        help="JSONL with accession and sequence per line; enables the sequence-level "
        "leakage check (default: skipped)",
    )
    parser.add_argument(
        "--recompute-auc",
        action="store_true",
        help="additionally recompute the AUC fields from the committed tables and "
        "models and cross-check metrics.json at 1e-9 (hard gate); this trains 200 "
        "models and takes minutes, so it is off by default",
    )
    parser.add_argument(
        "--strict-paper",
        action="store_true",
        help="promote a paper-layer miss beyond the 0.01 tolerance to exit 1",
    )
    parser.add_argument(
        "--strict-leakage",
        action="store_true",
        help="promote the train/test sequence-identity report row to exit 1",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        metavar="PATH",
        help="report path (default: <repo root>/validation_report.json)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="print the JSON report to stdout instead of the text table",
    )
    return parser.parse_args(argv)


def _load_sequences(path: Path) -> dict[str, str] | None:
    """Read an ``accession``/``sequence`` JSONL into a mapping, or ``None``.

    The committed tables carry no sequences and the cache is gitignored, so a
    missing file is not an error: the sequence-level check then SKIPs rather
    than passing (Ruling 2).
    """
    if not path.is_file():
        print(
            f"note: sequence cache {path} not found; the sequence-level leakage "
            "check will SKIP",
            file=sys.stderr,
        )
        return None
    sequences: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            record = json.loads(stripped)
            accession = record.get("accession")
            sequence = record.get("sequence")
            if accession is None or sequence is None:
                continue
            sequences[str(accession)] = str(sequence)
    return sequences


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    sequences = _load_sequences(args.sequences) if args.sequences is not None else None
    try:
        results = run_all(
            root=args.root,
            models_dir=args.models_dir,
            metrics_path=args.metrics,
            sequences=sequences,
            strict_paper=args.strict_paper,
            recompute_auc=args.recompute_auc,
        )
    except (ValidationError, TrainingDataSchemaError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    report = results_to_report(results)
    output = args.output if args.output is not None else REPO_ROOT / "validation_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if args.json:
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(format_report_table(results))
    return exit_code_for(
        results, strict_paper=args.strict_paper, strict_leakage=args.strict_leakage
    )


if __name__ == "__main__":
    raise SystemExit(main())
