
# -*- coding: utf-8 -*-
"""Stage 02: curation.

Every source entry is given an outcome. The scope decision is recorded
separately from the per-structure checks, being a decision about a family
rather than about a structure.
"""
from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

from .common import (CHEMISTRY_ERRORS, RunError, canonical, fragment_count,
                     has_metal, heavy_atoms, net_charge, parse,
                     radical_electrons)

#: The charge an ion of each type is required to carry. The scope of the method
#: is connected monovalent organic ions, so these are the only two values.
REQUIRED_CHARGE: Dict[str, int] = {"cation": 1, "anion": -1}

#: Families excluded by a decision about the family, taken before the run.
#: Neither is an organic ion of the kind the method edits: one is inorganic
#: anions, the other is coordination complexes of a metal centre.
OUT_OF_SCOPE_FAMILIES: Dict[str, str] = {
    "inorganics": "inorganic anion, not an organic ion",
    "metal complexes": "metal coordination complex, not an organic ion",
}

#: The order per-structure checks are applied in. Stated so the audit is
#: readable: an entry is reported against the first check it fails.
CHECK_ORDER: Sequence[str] = (
    "unreadable_structure",
    "disconnected",
    "radical",
    "metal_atom",
    "unknown_ion_type",
    "declared_charge_not_required",
    "computed_charge_not_required",
    "charge_not_on_a_single_atom",
)

def family_id(ion_type: str, family_raw: str) -> str:
    """A family is a (type, label) pair. Two types may share a label."""
    return "%s__%s" % (ion_type.strip().lower(),
                       family_raw.strip().lower().replace(" ", "_"))


def _structure_check(entry: Dict[str, Any]) -> Tuple[str, str, int]:
    """The first per-structure check the entry fails, with its canonical form.

    Returns (check, canonical_smiles, heavy_atoms). An empty check means the
    entry passed every per-structure test.
    """
    ion_type = entry["ion_type_raw"].strip().lower()
    try:
        mol = parse(entry["smiles_raw"])
    except CHEMISTRY_ERRORS:
        return "unreadable_structure", "", 0

    if fragment_count(mol) != 1:
        return "disconnected", "", heavy_atoms(mol)
    if radical_electrons(mol) != 0:
        return "radical", "", heavy_atoms(mol)
    if has_metal(mol):
        return "metal_atom", "", heavy_atoms(mol)
    if ion_type not in REQUIRED_CHARGE:
        return "unknown_ion_type", "", heavy_atoms(mol)

    required = REQUIRED_CHARGE[ion_type]
    declared = entry["declared_charge_raw"].strip()
    try:
        declared_value = int(float(declared))
    except (TypeError, ValueError):
        return "declared_charge_not_required", "", heavy_atoms(mol)
    if declared_value != required:
        return "declared_charge_not_required", "", heavy_atoms(mol)
    if net_charge(mol) != required:
        return "computed_charge_not_required", "", heavy_atoms(mol)
    # the charge has to be one atom's. A nitro, nitrate ester or azide group
    # carries a positive and a negative atom that cancel; the method holds one
    # charged atom and requires every other atom of a product to be neutral, so
    # an ion written this way is outside its scope.
    if sum(1 for a in mol.GetAtoms()
           if int(a.GetFormalCharge()) != 0) != 1:
        return "charge_not_on_a_single_atom", "", heavy_atoms(mol)

    return "", canonical(mol), heavy_atoms(mol)


def _assert_check_order() -> None:
    """CHECK_ORDER is the order the audit reports against, so it has to be the
    order the code applies. Any name returned has to appear in it."""
    import inspect
    source = inspect.getsource(_structure_check)
    named = [c for c in CHECK_ORDER if '"%s"' % c in source]
    missing = [c for c in CHECK_ORDER if c not in named]
    if missing:
        raise RunError("CHECK_ORDER names checks the code does not apply: %s"
                       % ", ".join(missing))


_assert_check_order()


def curate(entries: Sequence[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]],
                                                       List[Dict[str, Any]]]:
    """Return (audit rows, retained ions).

    One audit row per source entry. A retained ion carries an inventory id
    assigned in audit order, so the identifier is a property of the source
    ordering and of nothing else.
    """
    audit: List[Dict[str, Any]] = []
    retained: List[Dict[str, Any]] = []
    seen: Dict[Tuple[str, str], int] = {}

    for entry in entries:
        ion_type = entry["ion_type_raw"].strip().lower()
        fid = family_id(ion_type, entry["family_raw"])
        label = entry["family_raw"].strip().lower()

        in_scope = label not in OUT_OF_SCOPE_FAMILIES
        scope_reason = "" if in_scope else OUT_OF_SCOPE_FAMILIES[label]

        check, canon, heavies = _structure_check(entry)

        if not in_scope:
            status = "excluded_out_of_scope"
        elif check:
            status = "excluded_" + check
        else:
            key = (fid, canon)
            if key in seen:
                status = "excluded_duplicate_in_family"
            else:
                status = "retained"

        row = dict(entry)
        row.update({
            "family_id": fid,
            "required_charge": REQUIRED_CHARGE.get(ion_type, ""),
            "in_scope": int(in_scope),
            "scope_reason": scope_reason,
            "structure_check": check or "ok",
            "status": status,
            "canonical_smiles": canon,
            "heavy_atoms": heavies,
        })
        audit.append(row)

        if status == "retained":
            inventory_id = len(retained) + 1
            seen[(fid, canon)] = inventory_id
            retained.append({
                "inventory_id": inventory_id,
                "raw_index": entry["raw_index"],
                "family_id": fid,
                "family_label": entry["family_raw"].strip(),
                "ion_type": ion_type,
                "required_charge": REQUIRED_CHARGE[ion_type],
                "canonical_smiles": canon,
                "heavy_atoms": heavies,
                "abbreviation": entry["abbreviation_raw"],
                "chemical_name": entry["chemical_name_raw"],
            })
    return audit, retained


def summary(audit: Sequence[Dict[str, Any]],
            retained: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Counts, every one computed from the audit rows."""
    from collections import Counter
    status = Counter(r["status"] for r in audit)
    families = sorted({r["family_id"] for r in retained})
    by_family = Counter(r["family_id"] for r in retained)
    scope_excluded = Counter(r["family_raw"].strip().lower()
                             for r in audit if not r["in_scope"])
    return {
        "source_entries": len(audit),
        "retained_ions": len(retained),
        "families_retained": len(families),
        "status_counts": dict(sorted(status.items())),
        "out_of_scope_families": dict(sorted(scope_excluded.items())),
        "members_per_family": dict(sorted(by_family.items())),
        "check_order": list(CHECK_ORDER),
    }
