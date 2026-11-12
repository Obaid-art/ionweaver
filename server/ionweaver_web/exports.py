"""Files the user takes away: CSV, SDF and plain SMILES."""
from __future__ import annotations

import csv
import io
from typing import Any, Dict, Iterable, List

from rdkit import Chem

from .descriptors import DESCRIPTORS

__all__ = ["COLUMNS", "products_csv", "products_sdf", "products_smi"]

CHECK_KEYS = ("connected_closed_shell", "required_charge", "core_atoms", "core_bonds", "charge_on_core")

COLUMNS: List[str] = [
    "family", "ion_type", "parent_smiles", "product_smiles", "operator", "edit", "description", "site",
    "similarity_to_parent", "heavy_atom_change", "new_to_inventory", "certified",
    *("check_" + k for k in CHECK_KEYS), "judge_agrees", "definition_id", "grammar_id", "core_pattern",
    "formula", *(d["key"] for d in DESCRIPTORS), "mapped_smiles", "round", "rule",
]


def _flat(p: Dict[str, Any]) -> Dict[str, Any]:
    d = p.get("descriptors", {}) or {}
    row = {
        "family": p.get("family_label", ""), "ion_type": "cation" if str(p.get("family_id", "")).startswith("cation") else "anion",
        "parent_smiles": p.get("parent_smiles", ""), "product_smiles": p.get("product_smiles", ""),
        "operator": p.get("operator", ""), "edit": p.get("edit", ""), "description": p.get("description", ""),
        "site": p.get("site", ""), "similarity_to_parent": p.get("similarity_to_parent", ""),
        "heavy_atom_change": p.get("heavy_atom_change", ""), "new_to_inventory": int(bool(p.get("new_to_inventory"))),
        "certified": int(bool(p.get("certified"))), "judge_agrees": int(bool(p.get("judge_agrees"))),
        "definition_id": p.get("definition_id", ""), "grammar_id": p.get("grammar_id", ""),
        "core_pattern": p.get("core_pattern", ""), "formula": d.get("formula", ""),
        "mapped_smiles": p.get("mapped_smiles", ""), "round": p.get("round", ""), "rule": p.get("rule", ""),
    }
    for k in CHECK_KEYS:
        row["check_" + k] = int(bool((p.get("checks") or {}).get(k)))
    for desc in DESCRIPTORS:
        row[desc["key"]] = d.get(desc["key"], "")
    return row


def _safe(value: Any) -> Any:
    """A text cell that a spreadsheet would read as a formula gets a leading apostrophe."""
    if isinstance(value, str) and value[:1] in ("=", "+", "-", "@", "\t", "\r"):
        return "'" + value
    return value


def products_csv(products: Iterable[Dict[str, Any]]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for p in products:
        w.writerow({k: _safe(v) for k, v in _flat(p).items()})
    return buf.getvalue()


def products_smi(products: Iterable[Dict[str, Any]]) -> str:
    seen, lines = set(), []
    for p in products:
        smi = p.get("product_smiles", "")
        if smi and smi not in seen:
            seen.add(smi)
            lines.append(f"{smi}\t{p.get('family_label', '')}|{p.get('operator', '')}")
    return "\n".join(lines) + "\n"


def products_sdf(products: Iterable[Dict[str, Any]]) -> str:
    out = io.StringIO()
    writer = Chem.SDWriter(out)
    try:
        for p in products:
            mol = Chem.MolFromSmiles(str(p.get("product_smiles", "")))
            if mol is None:
                continue
            mol.SetProp("_Name", str(p.get("product_smiles", "")))
            for key, value in _flat(p).items():
                if value != "" and value is not None:
                    mol.SetProp(key, str(value))
            writer.write(mol)
    finally:
        writer.close()
    return out.getvalue()
