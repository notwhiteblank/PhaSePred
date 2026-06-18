#!/usr/bin/env python
"""Run PScore in parallel across N chunks of a sequence table."""
from __future__ import annotations

import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
PSCORE_DIR = REPO_ROOT / "Tools" / "PScore" / "SourceCodeS2"
PSCORE_SCRIPT = PSCORE_DIR / "elife_phase_separation_predictor.py"


def run_pscore_chunk(args: tuple[int, list[dict]]) -> tuple[int, str]:
    chunk_idx, records = args
    fasta_text = "".join(
        f">{r['UniprotEntry']}\n{str(r['sequence']).upper()}\n" for r in records
    )
    with tempfile.TemporaryDirectory(prefix=f"pscore_chunk{chunk_idx}_") as tmp:
        workdir = Path(tmp)
        fasta = workdir / "input.fasta"
        output = workdir / "pscore.tsv"
        fasta.write_text(fasta_text, encoding="utf-8")
        subprocess.run(
            [
                "uv", "run", "python", str(PSCORE_SCRIPT),
                str(fasta), "-output", str(output), "-overwrite", "-mute",
            ],
            cwd=PSCORE_DIR,
            check=True,
            capture_output=True,
            text=True,
        )
        return chunk_idx, output.read_text(encoding="utf-8", errors="replace")


def parse_pscore(text: str) -> list[dict[str, str | float]]:
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].startswith("PScore"):
            accession = parts[-1].lstrip(">")
            rows.append({"UniprotEntry": accession, "PScore": float(parts[1])})
    return rows


def main() -> None:
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} SEQUENCES_CSV OUTPUT_CSV [N_JOBS]", file=sys.stderr)
        sys.exit(2)

    seq_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    df = pd.read_csv(seq_path)
    records = (
        df[["UniprotEntry", "sequence"]]
        .dropna()
        .assign(UniprotEntry=lambda d: d["UniprotEntry"].astype(str).str.strip())
        .drop_duplicates("UniprotEntry")
        .to_dict(orient="records")
    )

    chunk_size = max(1, len(records) // n_jobs)
    chunks = []
    for i in range(0, len(records), chunk_size):
        chunks.append((len(chunks), records[i : i + chunk_size]))

    print(f"Running PScore over {len(records)} sequences in {len(chunks)} chunks...")
    all_rows = []
    with Pool(processes=min(len(chunks), n_jobs)) as pool:
        for chunk_idx, output_text in pool.imap_unordered(run_pscore_chunk, chunks):
            rows = parse_pscore(output_text)
            all_rows.extend(rows)
            print(f"  Chunk {chunk_idx + 1}/{len(chunks)}: {len(rows)} scores")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(all_rows)
    result.to_csv(output_path, index=False)
    print(f"Wrote {len(result)} PScore rows to {output_path}")


if __name__ == "__main__":
    main()
