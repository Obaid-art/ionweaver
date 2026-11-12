
# -*- coding: utf-8 -*-
"""Stage 07: sequential editing, twelve rounds from one ion per core.

The root for each core is the ion whose body size is closest to the median for
that core.

Three selection rules share one generator and one set of checks, and differ
only in which of the certified products is carried forward:

    random        any new valid product, with no distance preference
    local         the product most similar to the starting ion
    exploration   the product least similar to the starting ion

Similarity is measured against the ion the walk started from, not against the
ion of the current round. No rule applies a threshold.

Every candidate examined in every round is recorded, not only the one selected.

After round one the molecule being edited is a product, and its core is taken
from the certification label rather than from a fresh search on the product.
"""
from __future__ import annotations

import hashlib
import statistics
from dataclasses import dataclass
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

from .build import CORE, certify, core_induced_bonds, label, realise
from .common import STUDY_SCORING, AtomEnv, BondEnv, Scoring, parse
from .operators import Context, context_from, definitions_for
from .single import order_for, similarity

ROUNDS = 12
RULES = ("random", "local", "exploration")

#: what each rule does with the similarity of a product to the starting ion
PREFERENCE = {"random": "none", "local": "most similar",
              "exploration": "least similar"}

ROOT_COLUMNS = ("family_id", "grammar_id", "core_pattern", "inventory_id",
                "canonical_smiles", "core_atoms", "body_heavy_atoms",
                "median_body_heavy_atoms", "distance_from_median",
                "members_of_this_core", "chosen_by")

STEP_COLUMNS = ("family_id", "grammar_id", "rule", "preference", "round",
                "held_smiles", "held_heavy_atoms", "candidates_examined",
                "candidates_certified", "selected_smiles",
                "selected_similarity_to_held", "selected_similarity_to_start",
                "selected_operator", "selected_description", "held_or_moved",
                "reason")

CANDIDATE_COLUMNS = ("family_id", "grammar_id", "rule", "preference", "round",
                     "held_smiles", "operator", "mode", "product_smiles",
                     "similarity_to_held", "similarity_to_start", "certified",
                     "selected", "description")

SUMMARY_COLUMNS = ("family_id", "grammar_id", "rule", "preference",
                   "rounds_requested", "rounds_moved", "rounds_held",
                   "final_smiles", "final_heavy_atoms",
                   "similarity_first_to_last", "distinct_products",
                   "total_candidates_examined")


def roots(starting: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One ion per core: the member whose body is closest to the median."""
    groups: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    for record in starting:
        groups.setdefault((record["family_id"], record["grammar_id"]),
                          []).append(record)
    out: List[Dict[str, Any]] = []
    for key in sorted(groups):
        members = sorted(groups[key], key=lambda r: r["canonical_smiles"])
        sizes = [r["body_heavy_atoms"] for r in members]
        middle = statistics.median(sizes)
        best = min(members, key=lambda r: (abs(r["body_heavy_atoms"] - middle),
                                           r["canonical_smiles"]))
        out.append({
            "family_id": key[0], "grammar_id": key[1],
            "core_pattern": best["core_pattern"],
            "inventory_id": best["inventory_id"],
            "canonical_smiles": best["canonical_smiles"],
            "body_heavy_atoms": best["body_heavy_atoms"],
            "median_body_heavy_atoms": middle,
            "distance_from_median": abs(best["body_heavy_atoms"] - middle),
            "members_of_this_core": len(members),
            "chosen_by": "closest body size to the median for this core, "
                         "then canonical SMILES"})
    return out


def pick_key(seed: int, family_id: str, grammar_id: str, rule: str,
             round_number: int, smiles: str) -> str:
    """The keyed order a round's candidates are considered in.

    Under the random rule this order is the choice. Under the other two it is
    only the tie-break, so that two products equally similar to the starting
    ion are still chosen reproducibly.
    """
    return hashlib.blake2b(
        ("%d|%s|%s|%s|%d|%s" % (seed, family_id, grammar_id, rule,
                                round_number, smiles)).encode("utf-8"),
        digest_size=16).hexdigest()


def walk(root: Dict[str, Any], rule: str,
         library: Sequence[Dict[str, Any]],
         fragments: Sequence[Dict[str, Any]],
         observed: FrozenSet[AtomEnv], elements: Sequence[int], seed: int,
         allowance: int, observed_bonds: FrozenSet[BondEnv],
         scoring: Scoring = STUDY_SCORING
         ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]],
                    Dict[str, Any]]:
    """Twelve rounds from one root, under one selection rule."""
    preference = PREFERENCE[rule]
    required = 1 if root["family_id"].startswith("cation") else -1
    held = parse(root["canonical_smiles"])
    start = parse(root["canonical_smiles"])
    held_core = sorted(int(i) for i in str(root.get("core_atoms", "")).split())         if root.get("core_atoms") else None
    if held_core is None:
        raise ValueError("a root needs its core atoms")

    steps: List[Dict[str, Any]] = []
    candidates: List[Dict[str, Any]] = []
    products: Set[str] = set()
    first = Chem.MolToSmiles(held)
    moved = 0

    for round_number in range(1, ROUNDS + 1):
        held_smiles = Chem.MolToSmiles(held)
        ion = context_from(held, held_core, root["inventory_id"],
                           root["family_id"], held_smiles, elements)
        kept, _tally, _refused = definitions_for(ion, library, fragments,
                                                 observed, observed_bonds)
        labelled = label(ion.mol, sorted(ion.core))
        expected = core_induced_bonds(labelled, frozenset(ion.core))

        chosen = order_for(seed + round_number, root["inventory_id"],
                           "%s|%s" % (rule, root["grammar_id"]),
                           kept, allowance)
        found: List[Tuple[Dict[str, Any], Chem.Mol, List[int]]] = []
        for definition in chosen:
            product = realise(labelled, definition)
            if not product.built:
                continue
            certificate = certify(product.mol, labelled, required, expected,
                                  len(ion.core))
            near_held = similarity(ion.mol, product.mol)
            near_start = similarity(start, product.mol)
            row = {
                "family_id": root["family_id"],
                "grammar_id": root["grammar_id"], "rule": rule,
                "preference": preference, "round": round_number,
                "held_smiles": held_smiles, "operator": definition.operator,
                "mode": definition.mode, "product_smiles": product.smiles,
                "similarity_to_held": near_held,
                "similarity_to_start": near_start,
                "certified": int(certificate.certified),
                "selected": 0,
                "description": definition.description}
            candidates.append(row)
            if certificate.certified and product.smiles != held_smiles:
                core_now = sorted(
                    int(a.GetIdx()) for a in product.mol.GetAtoms()
                    if a.HasProp(CORE) and int(a.GetIntProp(CORE)) == 1)
                found.append((row, product.mol, core_now))

        examined = sum(1 for r in candidates if r["round"] == round_number
                       and r["held_smiles"] == held_smiles)
        certified = sum(1 for r in candidates if r["round"] == round_number
                        and r["held_smiles"] == held_smiles
                        and r["certified"] == 1)

        step = {"family_id": root["family_id"],
                "grammar_id": root["grammar_id"], "rule": rule,
                "preference": preference, "round": round_number,
                "held_smiles": held_smiles,
                "held_heavy_atoms": ion.mol.GetNumAtoms(),
                "candidates_examined": examined,
                "candidates_certified": certified}

        if not found:
            step.update({"selected_smiles": "",
                         "selected_similarity_to_held": "",
                         "selected_similarity_to_start": "",
                         "selected_operator": "",
                         "selected_description": "", "held_or_moved": "held",
                         "reason": "no certified product was offered"})
            steps.append(step)
            continue

        # every certified product is eligible under every rule. What differs
        # is only the order they are preferred in, and the keyed hash is the
        # tie-break so the choice is reproducible either way.
        def ordering(item):
            key = pick_key(seed, root["family_id"], root["grammar_id"], rule,
                           round_number, item[0]["product_smiles"])
            if preference == "most similar":
                return (-item[0]["similarity_to_start"], key)
            if preference == "least similar":
                return (item[0]["similarity_to_start"], key)
            return (0.0, key)

        found.sort(key=ordering)
        row, product_mol, core_now = found[0]
        row["selected"] = 1
        step.update({"selected_smiles": row["product_smiles"],
                     "selected_similarity_to_held": row["similarity_to_held"],
                     "selected_similarity_to_start":
                         row["similarity_to_start"],
                     "selected_operator": row["operator"],
                     "selected_description": row["description"],
                     "held_or_moved": "moved", "reason": ""})
        steps.append(step)
        products.add(row["product_smiles"])
        moved += 1
        held = product_mol
        held_core = core_now

    last = Chem.MolToSmiles(held)
    summary = {
        "family_id": root["family_id"], "grammar_id": root["grammar_id"],
        "rule": rule, "preference": preference, "rounds_requested": ROUNDS,
        "rounds_moved": moved, "rounds_held": ROUNDS - moved,
        "final_smiles": last, "final_heavy_atoms": held.GetNumAtoms(),
        "similarity_first_to_last": similarity(parse(first), parse(last)),
        "distinct_products": len(products),
        "total_candidates_examined": len(candidates)}
    return steps, candidates, summary
