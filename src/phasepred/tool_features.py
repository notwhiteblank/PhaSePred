"""Per-tool feature decomposition (the ``tool-features`` CLI surface).

Each PhaSePred feature is produced by an upstream tool (or a native/table
computation). This module computes them *individually* so a caller can
decompose a protein's PS-related profile tool-by-tool, matching the web
``tools.<x>`` payloads:

* summary value per tool (the ``single``/``NLLR``/``FCR`` column),
* optional per-residue arrays (``--residue-json``), for the tools whose local
  contracts expose them.

Tool name vocabulary is aligned with the web API (PLAN §6-S7): ``espritz,
seg, pscore, plaac, catgranule, deepcoil, deepphase, phos, hydropathy,
charge``. A single tool failing degrades to ``NaN`` for that tool only — the
whole run never aborts (same degrade semantics as ``predict``).
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import pandas as pd

from phasepred import tools
from phasepred.features import CHARGED_RESIDUES, compute_native_features

# Tool vocabulary (web API order).
TOOL_NAMES = [
    "espritz",
    "seg",
    "pscore",
    "plaac",
    "catgranule",
    "deepcoil",
    "deepphase",
    "phos",
    "hydropathy",
    "charge",
]

# Tools whose local contract exposes per-residue arrays for --residue-json.
_RESIDUE_TOOLS = {"espritz", "seg", "catgranule", "deepcoil", "charge"}


def _catgranule(records: list[dict[str, str]]) -> tuple[dict[str, float], dict[str, dict]]:
    from catgranule import score_sequence

    singles: dict[str, float] = {}
    residues: dict[str, dict] = {}
    for r in records:
        seq = str(r["sequence"]).upper()
        try:
            scored = score_sequence(seq)
            singles[r["accession"]] = float(scored["single"])
            residues[r["accession"]] = {"residue": [float(x) for x in scored["residue"]]}
        except Exception:
            singles[r["accession"]] = float("nan")
            residues[r["accession"]] = {"residue": []}
    return singles, residues


def _charge_residues(records: list[dict[str, str]]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in records:
        seq = str(r["sequence"]).upper()
        label = [1.0 if c in CHARGED_RESIDUES else 0.0 for c in seq]
        pos_pos = [i for i, c in enumerate(seq) if c in {"K", "R"}]
        pos_neg = [i for i, c in enumerate(seq) if c in {"D", "E"}]
        out[r["accession"]] = {"label": label, "POS_pos": pos_pos, "POS_neg": pos_neg}
    return out


def compute_tool_features(
    records: list[dict[str, str]],
    tool_names: list[str],
    *,
    residue: bool = False,
) -> tuple[pd.DataFrame, dict[str, dict[str, dict]]]:
    """Compute per-tool summary values (and optional residue arrays).

    Returns ``(summary_df, residue_map)``. ``summary_df`` has ``UniprotEntry``
    plus one column per requested tool (NaN on degrade). ``residue_map`` maps
    ``accession -> tool -> field -> list`` (empty when ``residue`` is False).
    """
    requested = [t for t in tool_names if t in TOOL_NAMES]
    accs = [r["accession"] for r in records]
    summary: dict[str, list[object]] = {t: [float("nan")] * len(accs) for t in requested}
    residue_map: dict[str, dict[str, dict]] = {a: {} for a in accs}

    # Native (hydropathy, charge) — pure Python, no tools.
    native: dict[str, dict[str, float]] = {a: {} for a in accs}
    for r in records:
        try:
            native[r["accession"]] = compute_native_features(str(r["sequence"]).upper())
        except Exception:
            native[r["accession"]] = {}

    if "hydropathy" in requested:
        for i, r in enumerate(records):
            summary["hydropathy"][i] = native[r["accession"]].get("Hydropathy", float("nan"))
    if "charge" in requested:
        for i, r in enumerate(records):
            summary["charge"][i] = native[r["accession"]].get("FCR", float("nan"))
        if residue:
            for a, ch_payload in _charge_residues(records).items():
                residue_map[a]["charge"] = ch_payload

    # catgranule (package) — pure Python.
    if "catgranule" in requested:
        singles, residues = _catgranule(records)
        for i, r in enumerate(records):
            summary["catgranule"][i] = singles[r["accession"]]
        if residue:
            for a, cg_payload in residues.items():
                residue_map[a]["catgranule"] = cg_payload

    # ESpritz (IDR) — perl tool; residue arrays derivable.
    if "espritz" in requested:
        esp_resid: dict[str, dict[str, object]] = (
            tools.espritz_residues(records) if residue else {}
        )
        if not residue:
            for i, acc in enumerate(accs):
                summary["espritz"][i] = tools.run_espritz(records).get(acc, float("nan"))
        else:
            for i, r in enumerate(records):
                esp_payload = esp_resid.get(r["accession"])
                if esp_payload is None:
                    summary["espritz"][i] = float("nan")
                    continue
                esp_label = cast(list[str], esp_payload["label"])
                summary["espritz"][i] = (
                    sum(1 for s in esp_label if s == "D") / len(esp_label)
                    if esp_label
                    else float("nan")
                )
                residue_map[r["accession"]]["espritz"] = {
                    "label": esp_label,
                    "residue": cast(list[float], esp_payload["residue"]),
                }

    # SEG (LCR) — native binary; residue mask derivable.
    if "seg" in requested:
        seg_resid: dict[str, dict[str, object]] = tools.seg_residues(records) if residue else {}
        if not residue:
            for i, acc in enumerate(accs):
                summary["seg"][i] = tools.run_seg(records).get(acc, float("nan"))
        else:
            for i, r in enumerate(records):
                seg_payload = seg_resid.get(r["accession"])
                if seg_payload is None:
                    summary["seg"][i] = float("nan")
                    continue
                seg_label = cast(list[float], seg_payload["label"])
                summary["seg"][i] = sum(seg_label) / len(seg_label) if seg_label else float("nan")
                residue_map[r["accession"]]["seg"] = {"label": seg_label}

    # DeepCoil (v2022 binary) — external env; residue arrays derivable.
    if "deepcoil" in requested:
        dc_resid: dict[str, dict[str, object]] = (
            tools.deepcoil_residues(records) if residue else {}
        )
        if not residue:
            for i, acc in enumerate(accs):
                summary["deepcoil"][i] = tools.run_deepcoil(records).get(acc, float("nan"))
        else:
            for i, r in enumerate(records):
                dc_payload = dc_resid.get(r["accession"])
                if dc_payload is None:
                    summary["deepcoil"][i] = float("nan")
                    continue
                raw_cc = cast(list[float], dc_payload["coiled-coil"])
                summary["deepcoil"][i] = (
                    float(1.0 if raw_cc and max(raw_cc) >= tools.DEEPCOIL_THRESHOLD else 0.0)
                    if raw_cc
                    else float("nan")
                )
                residue_map[r["accession"]]["deepcoil"] = {
                    "coiled-coil": raw_cc,
                    "sharpen": cast(list[float], dc_payload["sharpen"]),
                }

    # PScore / PLAAC — single only (residue not exposed by local contract).
    if "pscore" in requested:
        vals = tools.run_pscore(records)
        for i, r in enumerate(records):
            summary["pscore"][i] = vals.get(r["accession"], float("nan"))
    if "plaac" in requested:
        vals = tools.run_plaac(records)
        for i, r in enumerate(records):
            summary["plaac"][i] = vals.get(r["accession"], float("nan"))

    # Lookups (DeepPhase, PhosphoSitePlus).
    if "deepphase" in requested:
        vals = tools.lookup_deepphase(accs)
        for i, acc in enumerate(accs):
            summary["deepphase"][i] = vals.get(acc, float("nan"))
    if "phos" in requested:
        vals = tools.compute_phos_freq(records)
        for i, r in enumerate(records):
            summary["phos"][i] = vals.get(r["accession"], float("nan"))

    frame = pd.DataFrame({"UniprotEntry": accs})
    for t in requested:
        frame[t] = summary[t]
    return frame, residue_map


def write_residue_json(residue_map: dict[str, dict[str, dict]], out_dir: Path) -> None:
    """Write one JSON file per accession: ``{tool: {field: [...]}}``."""
    import json

    out_dir.mkdir(parents=True, exist_ok=True)
    for acc, tools_map in residue_map.items():
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in acc)
        path = out_dir / f"{safe}.json"
        path.write_text(json.dumps(tools_map, indent=2), encoding="utf-8")
