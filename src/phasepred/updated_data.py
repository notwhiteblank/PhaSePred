"""Build training label tables and splits from current (updated) databases.

PhaSepDB 3.0 provides PS-self / PS-other classification.
PhaSePro and LLPSDB2 positives are held out as external validation sets.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths relative to repo root
# ---------------------------------------------------------------------------

PHASEPDB3_UNIPROT = "data/raw/external/phasepdb3/uniprot.jsonl"
PHASEPRO_JSON = "data/raw/external/phasepro/download_full.json"
LLPSDB2_DIR = "data/raw/external/llpsdb2/extracted"

# LLPSDB2 subdirectory names (casing from actual filesystem)
LLPSDB2_POS = "Phase_separation_unambiguous/Phase_separation_Unambiguous"
LLPSDB2_NEG = "No_phase_separation_unambiguous/No_phase_separation_Unambiguous"


# ---------------------------------------------------------------------------
# Dataclass for a unified label record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LabelRecord:
    accession: str
    organism: str
    is_human: bool
    is_ps_self: bool
    is_ps_other: bool
    is_any_ps: bool
    is_ps_negative: bool
    source: str  # "phasepdb3", "phasepro", "llpsdb2"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_phasepdb3_labels(repo_root: str | Path) -> pd.DataFrame:
    """Load PhaSepDB 3.0 labels from uniprot.jsonl.

    Returns DataFrame with columns:
        UniprotEntry, organism, is_human, is_ps_self, is_ps_other,
        aggregated_class
    """
    path = Path(repo_root) / PHASEPDB3_UNIPROT
    rows: list[dict[str, object]] = []
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            acc = str(rec["uniprot_id"]).strip()
            organism = str(rec.get("organism", "")).strip()
            agg_class = str(rec.get("aggregated_class", "")).strip()

            rows.append(
                {
                    "UniprotEntry": acc,
                    "organism": organism,
                    "is_human": organism == "Homo sapiens",
                    "is_ps_self": "PS-self" in agg_class,
                    "is_ps_other": "PS-other" in agg_class,
                    "aggregated_class": agg_class,
                }
            )

    df = pd.DataFrame(rows)
    df["source"] = "phasepdb3"
    return df


def load_phasepro_labels(repo_root: str | Path) -> pd.DataFrame:
    """Load PhaSePro labels. All 121 entries are confirmed PS positives.

    Returns DataFrame with columns:
        UniprotEntry, organism, is_human, is_any_ps
    """
    path = Path(repo_root) / PHASEPRO_JSON
    with open(path) as fh:
        data = json.load(fh)

    rows: list[dict[str, object]] = []
    for acc, rec in data.items():
        organism = str(rec.get("organism", "")).strip()
        rows.append(
            {
                "UniprotEntry": str(acc).strip(),
                "organism": organism,
                "is_human": "Homo sapiens" in organism,
                "is_any_ps": True,
            }
        )

    df = pd.DataFrame(rows)
    df["source"] = "phasepro"
    return df


def load_llpsdb2_labels(repo_root: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load LLPSDB2 positive and negative label sets.

    Returns (positives_df, negatives_df), each with columns:
        UniprotEntry, organism, is_human
    """
    base = Path(repo_root) / LLPSDB2_DIR

    pos_path = base / LLPSDB2_POS / "protein.xls"
    neg_path = base / LLPSDB2_NEG / "protein.xls"

    pos_df = _read_llpsdb2_protein(pos_path)
    pos_df["is_any_ps"] = True
    pos_df["is_ps_negative"] = False
    pos_df["source"] = "llpsdb2"

    neg_df = _read_llpsdb2_protein(neg_path)
    neg_df["is_any_ps"] = False
    neg_df["is_ps_negative"] = True
    neg_df["source"] = "llpsdb2"

    return pos_df, neg_df


def _read_llpsdb2_protein(path: Path) -> pd.DataFrame:
    df = pd.read_excel(str(path))
    col = "Uniprot ID" if "Uniprot ID" in df.columns else "UniprotID"
    species_col = "Species"

    records = []
    for _, row in df.iterrows():
        acc = str(row.get(col, "")).strip()
        if not acc or acc == "nan":
            continue
        organism = str(row.get(species_col, "")).strip()
        records.append(
            {
                "UniprotEntry": acc,
                "organism": organism,
                "is_human": "Homo sapiens" in organism or organism == "Human",
            }
        )
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Label table builder
# ---------------------------------------------------------------------------

def build_label_table(repo_root: str | Path) -> pd.DataFrame:
    """Merge all label sources into a unified table.

    PhaSepDB3 labels take precedence on conflicts.
    LLPSDB2 confirmed negatives override when PhaSepDB3 has no opinion.
    """
    phasepdb = load_phasepdb3_labels(repo_root)
    phasepro = load_phasepro_labels(repo_root)
    llps_pos, llps_neg = load_llpsdb2_labels(repo_root)

    # Start with PhaSepDB3 as base (authoritative for self/other)
    label_table = phasepdb.copy()

    # Merge PhaSePro: add is_any_ps flag for held-out validation
    # If already in PhaSepDB3, PhaSepDB3's self/other labels stay
    pp_only = phasepro[~phasepro["UniprotEntry"].isin(label_table["UniprotEntry"])]
    if len(pp_only) > 0:
        pp_extra = pp_only.assign(
            is_ps_self=False, is_ps_other=False, aggregated_class=""
        )
        label_table = pd.concat([label_table, pp_extra], ignore_index=True)

    # Merge LLPSDB2 positives (held-out validation)
    llps_pos_only = llps_pos[~llps_pos["UniprotEntry"].isin(label_table["UniprotEntry"])]
    if len(llps_pos_only) > 0:
        llps_extra = llps_pos_only.assign(
            is_ps_self=False, is_ps_other=False, aggregated_class=""
        )
        label_table = pd.concat([label_table, llps_extra], ignore_index=True)

    # Merge LLPSDB2 negatives
    llps_neg_only = llps_neg[~llps_neg["UniprotEntry"].isin(label_table["UniprotEntry"])]
    if len(llps_neg_only) > 0:
        llps_extra = llps_neg_only.assign(
            is_ps_self=False, is_ps_other=False, aggregated_class=""
        )
        label_table = pd.concat([label_table, llps_extra], ignore_index=True)

    # Fill NaN bools
    for col in ["is_ps_self", "is_ps_other", "is_any_ps", "is_ps_negative"]:
        if col in label_table.columns:
            label_table[col] = label_table[col].fillna(False).astype(bool)
    for col in ["is_human"]:
        if col in label_table.columns:
            label_table[col] = label_table[col].fillna(False).astype(bool)

    return label_table


# ---------------------------------------------------------------------------
# Task split builder
# ---------------------------------------------------------------------------

def build_task_split(
    label_table: pd.DataFrame,
    recomputed_features_path: str | Path,
    task: str,
    test_size: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build train/test split for a single task from the label table.

    Args:
        label_table: unified label table from build_label_table()
        recomputed_features_path: path to recomputed_features_full.csv (for
            background negative pool)
        task: one of "SaPS", "PdPS", "hSaPS", "hPdPS"
        test_size: fraction for test (default 0.2)
        random_state: random seed

    Returns (train_df, test_df) each with columns:
        UniprotEntry, task, split, label
    """
    human_only = task.startswith("h")
    if task in ("SaPS", "hSaPS"):
        positive_col = "is_ps_self"
    else:
        positive_col = "is_ps_other"

    # --- Positives ---
    subset = label_table[label_table["is_human"] == True] if human_only else label_table
    positives = subset[subset[positive_col] == True]["UniprotEntry"].unique().tolist()

    # --- Negative pool ---
    # 1. PhaSepDB3 proteins that are not positive for this task
    phasepdb_neg = subset[
        (subset["source"] == "phasepdb3") & (subset[positive_col] == False)
    ]["UniprotEntry"].unique().tolist()

    # 2. LLPSDB2 confirmed negatives
    llps_neg = subset[subset["is_ps_negative"] == True]["UniprotEntry"].unique().tolist()

    # 3. Background from recomputed features (proteins not in any positive set)
    background = _load_background_negatives(
        recomputed_features_path, label_table, human_only,
    )

    # 4. Filter out any positives from negative pool
    pos_set = set(positives)
    neg_pool_all = phasepdb_neg + llps_neg + background
    neg_pool = [uid for uid in neg_pool_all if uid not in pos_set]

    # Deduplicate
    neg_pool = list(dict.fromkeys(neg_pool))

    # --- Split ---
    pos_train, pos_test = _stratified_split_list(positives, test_size, random_state)
    neg_train, neg_test = _stratified_split_list(neg_pool, test_size, random_state)

    train_rows = _make_split_rows(pos_train, neg_train, task, "train")
    test_rows = _make_split_rows(pos_test, neg_test, task, "test")

    return pd.DataFrame(train_rows), pd.DataFrame(test_rows)


def _load_background_negatives(
    recomputed_features_path: str | Path,
    label_table: pd.DataFrame,
    human_only: bool,
) -> list[str]:
    """Get background negatives from recomputed features CSV."""
    all_positive_ids = set()
    for col in ["is_ps_self", "is_ps_other"]:
        if col in label_table.columns:
            all_positive_ids.update(
                label_table[label_table[col] == True]["UniprotEntry"].tolist()
            )
    # Also exclude LLPSDB2 confirmed positives
    if "is_any_ps" in label_table.columns:
        all_positive_ids.update(
            label_table[label_table["is_any_ps"] == True]["UniprotEntry"].tolist()
        )

    background: list[str] = []
    with open(recomputed_features_path) as fh:
        header = fh.readline().strip().split(",")
        uniprot_idx = header.index("UniprotEntry")
        for line in fh:
            parts = line.strip().split(",")
            uid = parts[uniprot_idx]
            if uid not in all_positive_ids:
                background.append(uid)

    # For human tasks, further restrict to human entries.
    # DeepPhase column is only populated for human Swiss-Prot entries
    # (Phos freq is always 0.0 for non-human, so it doesn't filter).
    if human_only:
        human_background = []
        with open(recomputed_features_path) as fh:
            header_line = fh.readline().strip().split(",")
            deepphase_idx = header_line.index("DeepPhase")
            uniprot_idx = header_line.index("UniprotEntry")
            for line in fh:
                parts = line.strip().split(",")
                uid = parts[uniprot_idx]
                if uid in background:
                    val = parts[deepphase_idx] if deepphase_idx < len(parts) else ""
                    if val and val != "" and val.lower() != "nan":
                        try:
                            float(val)
                            human_background.append(uid)
                        except ValueError:
                            pass
            return human_background

    return background

    return background


def _stratified_split_list(
    items: list[str], test_size: float, random_state: int
) -> tuple[list[str], list[str]]:
    """Simple stratified-like split preserving order within each class."""
    rng = np.random.RandomState(random_state)
    indices = rng.permutation(len(items))
    n_test = max(1, int(len(items) * test_size))
    test_idx = set(indices[:n_test])
    train_idx = set(indices[n_test:])

    train = [items[i] for i in range(len(items)) if i in train_idx]
    test = [items[i] for i in range(len(items)) if i in test_idx]
    return train, test


def _make_split_rows(
    positives: list[str],
    negatives: list[str],
    task: str,
    split: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for uid in positives:
        rows.append({"UniprotEntry": uid, "task": task, "split": split, "label": 1})
    for uid in negatives:
        rows.append({"UniprotEntry": uid, "task": task, "split": split, "label": 0})
    return rows
