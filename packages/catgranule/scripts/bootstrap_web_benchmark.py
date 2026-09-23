"""S4 benchmark data bootstrap: parse the PhaSePred web-archive bulk dumps.

Reads the per-species ``<species>_reviewed.json`` bulk dumps
(``data/raw/external/phasepred_web/<file>.json``) and writes a compact
JSONL: one line per protein with ``accession``, ``entry_name``,
``organism``, ``sequence`` (the 2022-02-11 snapshot string used by the
archive), ``cg_single`` and ``cg_residue`` (the catGRANULE single value and
the comma-separated residue string parsed to a list).

The gate holdout split (SPEC 4.2) is also produced here
(``holdout_<species>.json``): a fixed-seed 80/20 split over the species'
protein list stored as accession lists.

Run: python packages/catgranule/scripts/bootstrap_web_benchmark.py
     [--species human mouse yeast] [--seed 42]

All outputs are gitignored (under data/raw/external/phasepred_web/).
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

OUT_DIR = Path("data/raw/external/phasepred_web")


def parse_species(species: str) -> tuple[Path, Path]:
    json_path = OUT_DIR / f"{species}_reviewed.json"
    compact_path = OUT_DIR / f"{species}_compact.jsonl"
    if not json_path.exists():
        raise SystemExit(f"Missing bulk dump: {json_path}")
    return json_path, compact_path


def load_species(species: str) -> dict[str, dict]:
    json_path = OUT_DIR / f"{species}_reviewed.json"
    print(f"[parse] loading {json_path} ...", flush=True)
    with json_path.open("rb") as fh:
        data = json.load(fh)
    return data


def write_compact(species: str, data: dict[str, dict], compact_path: Path) -> dict:
    n_cg = 0
    n_residue = 0
    residue_total = 0
    sequence_chars = 0
    errors = []
    with compact_path.open("w", encoding="utf-8") as out:
        for accession, rec in data.items():
            sequence = rec.get("Sequence", "")
            cat = rec.get("catGRANULE") or {}
            single = cat.get("single")
            residue_str = cat.get("residue")
            residue = None
            if isinstance(residue_str, str):
                try:
                    residue = [float(x) for x in residue_str.split(",") if x != ""]
                    n_residue += 1
                    residue_total += len(residue)
                except ValueError:
                    errors.append((accession, "residue parse"))
                    residue = None
            if single is not None:
                n_cg += 1
            sequence_chars += len(sequence)
            row = {
                "accession": accession,
                "entry_name": rec.get("Entry name"),
                "organism": rec.get("Organism"),
                "sequence": sequence,
                "cg_single": single,
                "cg_residue": residue,
            }
            out.write(json.dumps(row) + "\n")
    print(
        f"[parse] {species}: {len(data)} proteins, catgranule single={n_cg}, "
        f"residue arrays={n_residue} (total {residue_total} values, "
        f"{residue_total / max(n_residue, 1):.0f} avg len), {sequence_chars} seq chars",
        flush=True,
    )
    return {
        "entries": len(data),
        "entries_with_catgranule_single": n_cg,
        "entries_with_catgranule_residue": n_residue,
        "total_residue_values": residue_total,
        "errors": errors[:20],
        "n_errors": len(errors),
    }





def finalize_manifest(parsed: dict[str, dict], holdsplits: dict[str, dict]) -> None:
    raw_manifest = json.loads((OUT_DIR / "manifest_raw.json").read_text())
    for item in raw_manifest["files"]:
        species = item["organism"]
        item["parsed"] = parsed[species]
        item["holdout"] = {
            "n_eval_pool": holdsplits[species]["n_eval_pool"],
            "n_holdout": holdsplits[species]["n_holdout"],
            "file": f"holdout_{species}.json",
        }
    out = {
        "collection": "phasepred-web-archive",
        "data_version": "2022-02-11",
        "docs": "docs/CATGRANULE_SPEC.md 4.2 / docs/FEATURES.md appendix B",
        "units": "flux dict keyed by UniProt accession; 'catGRANULE'; residue csv strings",
        "files": raw_manifest["files"],
        "holdout_rule": "SPEC 4.2 shuffled seed 42; first 20% = gate holdout",
        "holdout_rule_seed": "42",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(out, indent=2), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=["human", "mouse", "yeast"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-parse", action="store_true")
    args = ap.parse_args()

    parsed: dict[str, dict] = {}
    holdsplits: dict[str, dict] = {}
    for species in args.species:
        json_path, compact_path = parse_species(species)
        if args.skip_parse and compact_path.exists():
            parsed[species] = {"skipped": True, "compact_file": compact_path.name}
            print(f"[parse] using existing {compact_path}", flush=True)
        else:
            data = load_species(species)
            parsed[species] = write_compact(species, data, compact_path)
            del data
        split_file = OUT_DIR / f"holdout_{species}.json"
        if split_file.exists() and not args.skip_parse:
            holdsplits[species] = json.loads(split_file.read_text())
            print(f"[holdout] {species}: using existing {split_file}", flush=True)
        else:
            accessions = []
            with compact_path.open() as fh:
                for line in fh:
                    accessions.append(json.loads(line)["accession"])
            holdsplits[species] = _split_from_accessions(species, accessions, args.seed)
    finalize_manifest(parsed, holdsplits)
    print("wrote manifest.json + compact files + holdouts", flush=True)
    return 0


def _split_from_accessions(species: str, accessions: list[str], seed: int) -> dict:
    rng = random.Random(seed)
    accs = sorted(accessions)
    rng.shuffle(accs)
    n = max(1, int(round(len(accs) * 0.2)))
    holdout = accs[:n]
    eval_pool = accs[n:]
    assert not (set(holdout) & set(eval_pool))
    out = {
        "species": species,
        "seed": seed,
        "n_proteins": len(accs),
        "n_holdout": len(holdout),
        "n_eval_pool": len(eval_pool),
        "fraction_holdout": 0.2,
        "holdout": holdout,
        "eval_pool": eval_pool,
    }
    (OUT_DIR / f"holdout_{species}.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(
        f"[holdout] {species}: {len(accs)} -> {len(eval_pool)} eval + "
        f"{len(holdout)} holdout (seed={seed})",
        flush=True,
    )
    return out


if __name__ == "__main__":
    sys.exit(main())