
# -*- coding: utf-8 -*-
"""Stage 03 — the layered discovery loop.

Each layer works on the members no accepted core has matched yet. The cores a
family accumulates are numbered G0, G1 and G2 in the order they were accepted,
so G0 is the motif the most members share and each later one holds what the
earlier ones did not reach.

Three candidate passes run per layer, on two ladders: the whole residual down a
long ladder, each charge-environment group down a short one, and every pair
once the residual is small.
"""
from __future__ import annotations

import collections
import itertools
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

from .common import parse
from .cores import (MIN_SUPPORT, PAIRWISE_RESIDUAL_LIMIT, candidate_thresholds,
                    charge_environment, charge_environment_thresholds,
                    core_query_key, evaluate_candidate, mcs_result)


def _raw_candidates(residual: Sequence[Dict[str, Any]],
                    mols: Dict[int, Chem.Mol]) -> List[Tuple[str, float, Chem.Mol]]:
    """The three passes, in order."""
    found: List[Tuple[str, float, Chem.Mol]] = []
    residual_mols = [mols[r["inventory_id"]] for r in residual]

    for threshold in candidate_thresholds(len(residual)):
        query = mcs_result(residual_mols, threshold)
        if query is not None:
            found.append(("family_residual", threshold, query))

    clusters: Dict[str, List[Dict[str, Any]]] = collections.defaultdict(list)
    for record in residual:
        clusters[charge_environment(mols[record["inventory_id"]])].append(record)
    for key in sorted(clusters):
        group = clusters[key]
        if len(group) < 2:
            continue
        group_mols = [mols[r["inventory_id"]] for r in group]
        for threshold in charge_environment_thresholds(len(group)):
            query = mcs_result(group_mols, threshold)
            if query is not None:
                found.append(("charge_environment", threshold, query))

    if len(residual) <= PAIRWISE_RESIDUAL_LIMIT:
        for a, b in itertools.combinations(residual, 2):
            query = mcs_result([mols[a["inventory_id"]],
                                mols[b["inventory_id"]]], 1.0)
            if query is not None:
                found.append(("pairwise", 1.0, query))
    return found


def discover(inventory: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Every family's cores, with the record of every decision taken."""
    families: Dict[str, List[Dict[str, Any]]] = collections.defaultdict(list)
    for record in inventory:
        families[record["family_id"]].append(record)

    models: List[Dict[str, Any]] = []
    assignments: List[Dict[str, Any]] = []
    diagnostics: List[Dict[str, Any]] = []
    decisions: List[Dict[str, Any]] = []
    uncovered: List[Dict[str, Any]] = []

    for family_id in sorted(families):
        members = sorted(families[family_id], key=lambda r: r["inventory_id"])
        mols = {r["inventory_id"]: parse(r["canonical_smiles"])
                for r in members}
        residual = list(members)
        layer = 0

        while len(residual) >= MIN_SUPPORT:
            seen: set = set()
            candidates: List[Tuple[Tuple, Any, str]] = []

            for source_pass, threshold, raw in _raw_candidates(residual, mols):
                key = core_query_key(raw)
                if key in seen:
                    continue
                seen.add(key)

                row = {"family_id": family_id, "layer": layer,
                       "pass": source_pass, "support_threshold": threshold,
                       "search_pattern": key, "pattern": "", "atoms": "",
                       "bonds": "", "charged_atoms": "", "support": 0,
                       "form": "", "copied_from": "",
                       "status": "rejected", "reason": ""}
                evaluated, declined = evaluate_candidate(raw, residual, mols)
                row["reason"] = declined
                if evaluated is not None:
                    query = evaluated.query
                    charged = sum(1 for a in query.GetAtoms()
                                  if int(a.GetFormalCharge()) != 0)
                    score = (len(evaluated.support), query.GetNumAtoms(),
                             query.GetNumBonds(), charged)
                    smarts = core_query_key(query)
                    row.update({"pattern": smarts,
                                "atoms": query.GetNumAtoms(),
                                "bonds": query.GetNumBonds(),
                                "charged_atoms": charged,
                                "support": len(evaluated.support),
                                "form": evaluated.form,
                                "copied_from": evaluated.copied_from,
                                "status": "candidate", "reason": ""})
                    candidates.append((score, evaluated, smarts))
                diagnostics.append(row)

            if not candidates:
                break

            candidates.sort(key=lambda c: (c[0], c[2]), reverse=True)
            score, chosen, smarts = candidates[0]
            grammar_id = "G%d" % sum(1 for m in models
                                     if m["family_id"] == family_id)

            def record_decision(point: str, rule: str, chosen_value: str,
                                trace: Dict[str, Any]) -> None:
                decisions.append({"family_id": family_id, "layer": layer,
                                  "grammar_id": grammar_id,
                                  "decision_point": point, "rule": rule,
                                  "chosen": chosen_value, **trace})

            runner_up = core_query_key(candidates[1][1].query) \
                if len(candidates) > 1 else ""
            tied = sum(1 for c in candidates if c[0] == score)
            record_decision(
                "4_which_candidate_becomes_the_core",
                "members, atoms, bonds, charged atoms, then pattern string",
                smarts,
                {"runner_up": runner_up,
                 "key_that_discriminated": ("pattern_string" if tied > 1
                                            else "members_atoms_bonds_charge"),
                 "options_tied_at_that_key": tied,
                 "options_considered": len(candidates)})
            record_decision(
                "2_which_member_the_pattern_is_copied_from",
                "fewest heavy atoms, then canonical SMILES",
                str(chosen.copied_from), chosen.source_trace)
            record_decision(
                "3_which_form_of_the_pattern_is_carried",
                "the form preserving the most support, then more atoms, "
                "then more bonds",
                chosen.form, chosen.form_trace)

            models.append({
                "family_id": family_id, "grammar_id": grammar_id,
                "layer": layer, "pattern": smarts,
                "atoms": chosen.query.GetNumAtoms(),
                "bonds": chosen.query.GetNumBonds(),
                "charged_atoms": score[3],
                "support": len(chosen.support),
                "support_fraction_of_family": round(
                    len(chosen.support) / len(members), 4),
                "copied_from_inventory_id": chosen.copied_from,
                "query_form": chosen.form})

            for inventory_id, match in sorted(chosen.support.items()):
                mol = mols[inventory_id]
                # the core is the match: the pattern that was found is the
                # pattern that is held, on every member of every family
                trace = chosen.traces[inventory_id]
                assignments.append({
                    "inventory_id": inventory_id, "family_id": family_id,
                    "grammar_id": grammar_id, "pattern": smarts,
                    "matched_atoms": " ".join(str(i) for i in match),
                    "usable_matches_available": trace["options_considered"],
                    "match_chosen_by_key": trace["key_that_discriminated"]})
                record_decision(
                    "1_which_atoms_of_a_member_are_its_core",
                    "bonds leaving, heteroatoms, aromatic atoms, index order",
                    "ion %d: %s" % (inventory_id, list(match)), trace)

            covered = set(chosen.support)
            residual = [r for r in residual
                        if r["inventory_id"] not in covered]
            layer += 1

        for record in residual:
            uncovered.append({
                "inventory_id": record["inventory_id"],
                "family_id": family_id,
                "canonical_smiles": record["canonical_smiles"],
                "reason": ("fewer than %d members remained" % MIN_SUPPORT
                           if len(residual) < MIN_SUPPORT
                           else "no candidate reached the floors or the "
                                "support minimum")})

    return {"models": models, "assignments": assignments,
            "diagnostics": diagnostics, "decisions": decisions,
            "uncovered": uncovered,
            "settings": {"min_support": MIN_SUPPORT,
                         "min_core_atoms": 3, "min_core_bonds": 2,
                         "pairwise_residual_limit": PAIRWISE_RESIDUAL_LIMIT,
                         "mcs_timeout": 0,
                         "bond_comparator": "CompareOrder"}}
