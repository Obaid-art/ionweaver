
# -*- coding: utf-8 -*-
"""The worked-example table.

Every term the paper defines is shown
working on one real ion, and the paper's explanations and figures are built from
this file and from nothing else.

The selection is deterministic and the rule is stated in every row. The
illustrating case is the smallest ion that exhibits it, breaking ties on
canonical SMILES. An example chosen by hand cannot be regenerated and is not
acceptable.

Everything here is read out of the archive rather than recomputed, for the same
reason `FACTS.json` is: an example in a figure and a row in the deposit cannot
then disagree.
"""
from __future__ import annotations

import csv
import os
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

COLUMNS = (
    "example_id", "illustrates", "kind", "inventory_id", "family_id",
    "core_pattern", "operator", "mode", "site", "atoms_or_bonds",
    "quantity_before", "quantity_after", "parent_smiles", "parent_selfies",
    "product_smiles", "product_selfies", "mapped_product_smiles",
    "selection_rule", "description")

#: The rule every row is chosen by, written into the row.
RULE = ("the smallest ion that exhibits the case, breaking ties on canonical "
        "SMILES")


def read(root: str, relative: str) -> List[Dict[str, str]]:
    path = os.path.join(root, relative)
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def heavy_atoms(smiles: str) -> str:
    """The heavy-atom count of a recorded SMILES.

    Derived from the file, not carried forward from the stage that wrote it,
    which is the same discipline `FACTS.json` follows.
    """
    if not smiles:
        return ""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    return str(mol.GetNumAtoms()) if mol is not None else ""


def _size(row: Dict[str, str]) -> Tuple[int, str]:
    """The key the selection rule sorts on."""
    parent = row.get("parent_smiles") or row.get("canonical_smiles") or ""
    return (len(parent), parent)


def smallest(rows: Sequence[Dict[str, str]],
             where: Callable[[Dict[str, str]], bool]
             ) -> Optional[Dict[str, str]]:
    found = [r for r in rows if where(r)]
    return min(found, key=_size) if found else None


def _from_examination(row: Dict[str, str], illustrates: str, kind: str,
                      before: Any = "", after: Any = "") -> Dict[str, Any]:
    return {
        "illustrates": illustrates, "kind": kind,
        "inventory_id": row["inventory_id"], "family_id": row["family_id"],
        "core_pattern": row["core_pattern"], "operator": row["operator"],
        "mode": row["mode"], "site": row["site"],
        "atoms_or_bonds": row["parameters"],
        "quantity_before": before, "quantity_after": after,
        "parent_smiles": row["parent_smiles"],
        "parent_selfies": row["parent_selfies"],
        "product_smiles": row["product_smiles"],
        "product_selfies": row["product_selfies"],
        "mapped_product_smiles": row["mapped_product_smiles"],
        "selection_rule": RULE, "description": row["description"]}


def build(root: str) -> List[Dict[str, Any]]:
    """Every row E3 asks for, chosen by the stated rule."""
    examinations = read(root, "05_single_edits/examination_records.csv")
    refusals = read(root, "05_single_edits/definition_refusals.csv")
    decisions = read(root, "03_cores/core_selection_decisions.csv")
    assignments = read(root, "03_cores/source_core_assignments.csv")
    perception = read(root, "09_validation/perception_check.csv")
    negatives = read(root, "09_validation/negative_tests.csv")
    out: List[Dict[str, Any]] = []

    certified = [r for r in examinations if r.get("certified") == "1"]

    # 1. each operator and each of its modes
    for operator, mode in sorted({(r["operator"], r["mode"])
                                  for r in certified}):
        row = smallest(certified, lambda r, o=operator, m=mode:
                       r["operator"] == o and r["mode"] == m)
        if row:
            out.append(_from_examination(
                row, "the %s operator, %s mode" % (operator, mode),
                "operator",
                before="%s heavy atoms in the parent"
                       % heavy_atoms(row["parent_smiles"]),
                after="%s heavy atoms in the product, a change of %s"
                      % (heavy_atoms(row["product_smiles"]),
                         row["heavy_atom_change"])))

    # 2. the definition-stage conditions, each shown declining something
    for condition, label in (
            ("ledger", "the valence ledger, declining an edit that has no "
                       "legal valence state"),
            ("atom_environment", "the atom-environment condition, declining a "
                                 "resulting environment no curated ion "
                                 "displays"),
            ("bond_environment", "the bond-environment condition, declining "
                                 "a bond between two kinds of atom no curated "
                                 "ion joins at that order"),
            ("neutral_and_non_aromatic", "the requirement that a constructed "
                                         "atom be neutral and non-aromatic"),
            ("terminal_halogen", "the terminal-halogen restriction"),
            ("above_lowest_valence_with_hydrogen",
             "the restriction on an atom above its lowest standard valence "
             "while it still carries hydrogen")):
        found = [r for r in refusals if r["condition"] == condition]
        if not found:
            out.append({
                "illustrates": label, "kind": "condition",
                "inventory_id": "", "family_id": "", "core_pattern": "",
                "operator": "", "mode": "", "site": "",
                "atoms_or_bonds": "", "quantity_before": "",
                "quantity_after": "0 definitions declined",
                "parent_smiles": "", "parent_selfies": "",
                "product_smiles": "", "product_selfies": "",
                "mapped_product_smiles": "", "selection_rule": RULE,
                "description": ("this condition declined no definition in the "
                                "run; an earlier condition reaches every case "
                                "it would have caught, and the shrinkage table "
                                "reports the zero")})
            continue
        pick = max(found, key=lambda r: int(r["definitions_refused"]))
        out.append({
            "illustrates": label, "kind": "condition",
            "inventory_id": pick["inventory_id"],
            "family_id": pick["family_id"], "core_pattern": "",
            "operator": pick["operator"], "mode": pick["mode"],
            "site": pick["example_site"],
            "atoms_or_bonds": pick["resulting_environment"],
            "quantity_before": "%s definitions reached this condition"
                               % pick["definitions_refused"],
            "quantity_after": "0 of them were created",
            "parent_smiles": "", "parent_selfies": "", "product_smiles": "",
            "product_selfies": "", "mapped_product_smiles": "",
            "selection_rule": "the row declining the most definitions under "
                              "this condition",
            "description": "%s. %s" % (pick["example_description"],
                                       pick["reason"])})

    # 3. a ring that arrived intact inside a fragment, matched to its source
    carried = [r for r in perception
               if r.get("source", "").startswith("the fragment")
               and r.get("agrees") == "1"]
    if carried:
        pick = min(carried, key=lambda r: (int(r["ring_size"]),
                                           r["inventory_id"]))
        row = smallest(certified,
                       lambda r, p=pick: r["examination_id"]
                       == p["examination_id"])
        base = _from_examination(
            row, "a ring that arrived intact inside the fragment an edit "
                 "introduced", "ring invariant",
            before="%s-ring of %s in the fragment"
                   % (pick["source_ring_size"], pick["source_bond_orders"]),
            after="%s-ring of %s in the product"
                  % (pick["ring_size"], pick["bond_orders"])) if row else None
        if base:
            base["description"] = (
                "the ring came from %s and is unchanged in size and in bond "
                "orders; no operator built it" % pick["source"])
            out.append(base)

    # 4. the five certification requirements, met and, where possible, failed
    for row in negatives:
        out.append({
            "illustrates": "certification requirement %s, %s"
                           % (row.get("requirement", ""),
                              "accepted" if row.get("should_certify") == "1"
                              else "refused"),
            "kind": "certification",
            "inventory_id": row.get("inventory_id", ""),
            "family_id": row.get("family_id", ""),
            "core_pattern": row.get("core_pattern", ""),
            "operator": "", "mode": "", "site": row.get("construction", ""),
            "atoms_or_bonds": row.get("atoms_or_bonds", ""),
            "quantity_before": row.get("expected", ""),
            "quantity_after": row.get("outcome", ""),
            "parent_smiles": row.get("parent_smiles", ""),
            "parent_selfies": "", "product_smiles": row.get("smiles", ""),
            "product_selfies": "", "mapped_product_smiles": "",
            "selection_rule": row.get("selection_rule", RULE),
            "description": row.get("note", "")})

    # 5. each core-discovery tie-break key, with the options that were tied
    for point in sorted({r["decision_point"] for r in decisions}):
        found = [r for r in decisions if r["decision_point"] == point]
        tied = [r for r in found
                if int(r.get("options_tied_at_that_key") or 0) > 1]
        pick = (min(tied, key=lambda r: (-int(r["options_tied_at_that_key"]),
                                         r["family_id"]))
                if tied else min(found, key=lambda r: r["family_id"]))
        out.append({
            "illustrates": "core-discovery decision point %s" % point,
            "kind": "core selection", "inventory_id": "",
            "family_id": pick["family_id"],
            "core_pattern": pick.get("chosen", ""), "operator": "",
            "mode": pick.get("grammar_id", ""), "site": pick.get("rule", ""),
            "atoms_or_bonds": pick.get("key_that_discriminated", ""),
            "quantity_before": "%s options considered"
                               % pick.get("options_considered", ""),
            "quantity_after": "%s still tied at the discriminating key"
                              % pick.get("options_tied_at_that_key", ""),
            "parent_smiles": "", "parent_selfies": "", "product_smiles": "",
            "product_selfies": "", "mapped_product_smiles": "",
            "selection_rule": "the decision with the most options tied at its "
                              "discriminating key, then the family name",
            "description": "the rule applied was: %s; the runner-up was %s"
                           % (pick.get("rule", ""),
                              pick.get("runner_up", "") or "none")})

    # 6. one accepted and one declined definition per operator
    for operator in sorted({r["operator"] for r in certified}):
        accepted = smallest(certified,
                            lambda r, o=operator: r["operator"] == o)
        if accepted:
            out.append(_from_examination(
                accepted, "an accepted definition of %s" % operator,
                "accepted definition",
                after="certified, similarity %s to the parent"
                      % accepted["similarity_to_parent"]))
        declined = [r for r in refusals if r["operator"] == operator]
        if declined:
            pick = max(declined, key=lambda r: int(r["definitions_refused"]))
            out.append({
                "illustrates": "a declined definition of %s" % operator,
                "kind": "declined definition",
                "inventory_id": pick["inventory_id"],
                "family_id": pick["family_id"], "core_pattern": "",
                "operator": operator, "mode": pick["mode"],
                "site": pick["example_site"],
                "atoms_or_bonds": pick["resulting_environment"],
                "quantity_before": "", "quantity_after": "never created",
                "parent_smiles": "", "parent_selfies": "",
                "product_smiles": "", "product_selfies": "",
                "mapped_product_smiles": "",
                "selection_rule": "the condition declining the most "
                                  "definitions of this operator",
                "description": "%s was declined by %s: %s"
                               % (pick["example_description"],
                                  pick["condition"], pick["reason"])})

    for number, row in enumerate(out, start=1):
        row["example_id"] = "E%03d" % number
    return out
