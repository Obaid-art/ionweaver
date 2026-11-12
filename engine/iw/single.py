
# -*- coding: utf-8 -*-
"""Stage 05 — single editing, under an allowance that every examination spends.

Up to the allowance's worth of candidates are examined per starting ion
per operator, and every examination spends one unit whatever it produces: a
refusal, a product identical to the parent, and a repeat each cost the same as a
new product. Where an operator offers fewer than the allowance, all are examined
and the number is recorded.

The order is a keyed shuffle and not an affine progression. The previous
order was a lattice, which can systematically over-represent library components
that happen to sit next to each other. Here each candidate is given a key from
the run seed, the ion, the operator and the candidate's own signature, and the
candidates are taken in key order. Because the key does not depend on the
allowance, a smaller allowance is a prefix of a larger one by construction --
and the run verifies it anyway, since a construction argument is not a test.

Keying on the signature rather than on the enumeration index means the order
survives a change to the order the enumerators yield in, so a rerun after an
unrelated edit to the code produces the same examinations.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from .common import (CHEMISTRY_ERRORS, STUDY_SCORING, Scoring, Verdict,
                     canonical, judge, parse)
from .build import (Certificate, Product, certify, core_induced_bonds, label,
                    realise)
from .perception import check as perceive
from .operators import Context, Definition

#: The allowances D3 asks for. The study reports the middle one; the other two
#: are the sensitivity check, and each must be a prefix of the next.
ALLOWANCES: Tuple[int, ...] = (25, 50, 100)

#: Similarity is Tanimoto over Morgan fingerprints at this radius and width.
#: Recorded in the settings file, because a similarity figure means nothing
#: without them.
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=FINGERPRINT_RADIUS,
                                                 fpSize=FINGERPRINT_BITS)

PT = Chem.GetPeriodicTable()

EXAMINATION_COLUMNS = (
    "examination_id", "inventory_id", "family_id", "core_pattern",
    "grammar_id", "operator", "mode", "order_index", "site", "definition_id",
    "parameters", "heavy_atom_change", "certified", "failing_requirement",
    "failing_reason", "parent_smiles", "parent_selfies", "product_smiles",
    "product_selfies", "mapped_product_smiles", "identical_to_parent",
    "repeat_for_this_operator", "repeat_for_this_ion",
    "similarity_to_parent", "judge_accepted", "judge_reason",
    "judge_agrees_with_certification", "hydrogens_written_explicitly",
    "description")

DEFINITION_SIGNATURE_COLUMNS = (
    "definition_id", "inventory_id", "family_id", "operator", "mode",
    "signature", "description")

LEDGER_COLUMNS = (
    "examination_id", "inventory_id", "family_id", "operator", "mode", "site",
    "atom", "role", "element", "formal_charge", "aromatic",
    "bond_orders_before", "hydrogens_before", "valence_spent_before",
    "bond_orders_after", "hydrogens_after", "valence_spent_after",
    "valence_spent_change", "total_valence_before", "total_valence_after",
    "answers_to", "environment_after")

REFUSAL_COLUMNS = ("examination_id", "inventory_id", "family_id", "operator",
                   "mode", "requirement", "reason", "product_smiles",
                   "description")


def fingerprint(mol: Chem.Mol):
    return _GEN.GetFingerprint(mol)


def similarity(parent: Chem.Mol, product: Chem.Mol) -> float:
    return round(DataStructs.TanimotoSimilarity(fingerprint(parent),
                                                fingerprint(product)), 6)


def selfies_of(smiles: str) -> str:
    """The SELFIES form, or an empty string where the encoder declines.

    Recorded rather than relied on. The string controls of D5 need it; nothing
    in the graph method does.
    """
    try:
        import selfies
        return selfies.encoder(smiles)
    except Exception:                                      # noqa: BLE001
        return ""


def shuffle_key(seed: int, inventory_id: int, operator: str,
                signature: str) -> str:
    """The key one candidate is taken in.

    A keyed hash, so the order is reproducible from the seed and independent of
    the order the enumerators happen to yield in, and so no arithmetic
    progression can line up with the library's own ordering.
    """
    digest = hashlib.blake2b(
        ("%d|%d|%s|%s" % (seed, inventory_id, operator, signature))
        .encode("utf-8"), digest_size=16)
    return digest.hexdigest()


def order_for(seed: int, inventory_id: int, operator: str,
              definitions: Sequence[Definition],
              allowance: Optional[int]) -> List[Definition]:
    """The candidates to examine, in the order they are examined."""
    keyed = sorted(definitions,
                   key=lambda d: (shuffle_key(seed, inventory_id, operator,
                                              d.signature_text),
                                  d.signature_text))
    return keyed if allowance is None else keyed[:allowance]


def prefix_check(seed: int, inventory_id: int, operator: str,
                 definitions: Sequence[Definition]) -> Tuple[bool, str]:
    """Is every smaller allowance a strict prefix of every larger one?"""
    orders = {n: [d.signature_text for d in
                  order_for(seed, inventory_id, operator, definitions, n)]
              for n in ALLOWANCES}
    for small, large in zip(ALLOWANCES, ALLOWANCES[1:]):
        shorter, longer = orders[small], orders[large]
        if shorter != longer[:len(shorter)]:
            for position, (a, b) in enumerate(zip(shorter, longer)):
                if a != b:
                    return False, ("allowance %d and allowance %d differ at "
                                   "position %d" % (small, large, position))
            return False, ("allowance %d is not a prefix of allowance %d"
                           % (small, large))
    return True, ""


def examine_ion(ion: Context, definitions: Sequence[Definition],
                grammar_id: str, core_pattern: str, required_charge: int,
                seed: int, scoring: Scoring = STUDY_SCORING,
                allowance: Optional[int] = None, first_id: int = 1
                ) -> Tuple[List[Dict[str, Any]], Dict[str, Any], List[str],
                           List[Dict[str, Any]], List[Dict[str, Any]],
                           List[Dict[str, Any]]]:
    """Examine one ion's candidates, operator by operator.

    Returns the examination rows, a per-ion-per-operator summary, any problem
    that should stop the run, and the ring and aromaticity rows the perception
    check of C2 produced.

    The perception check runs here, on the product while it is still in hand
    and still carrying the labels saying where each of its atoms came from.
    Run as a later pass it would have to rebuild every product to recover
    those labels, and a check that needs the thing it checks to be rebuilt is
    a check on a rebuild.
    """
    if allowance is None:
        allowance = ALLOWANCES[1] if scoring.allowance is None \
            else scoring.allowance

    labelled = label(ion.mol, sorted(ion.core))
    expected = core_induced_bonds(labelled, frozenset(ion.core))
    query = Chem.MolFromSmarts(core_pattern) if core_pattern else None
    parent_smiles = ion.canonical_smiles
    parent_selfies = selfies_of(parent_smiles)
    parent_mol = parse(parent_smiles)

    by_operator: Dict[str, List[Definition]] = {}
    for definition in definitions:
        by_operator.setdefault(definition.operator, []).append(definition)

    rows: List[Dict[str, Any]] = []
    summary: Dict[str, Any] = {}
    problems: List[str] = []
    ring_rows: List[Dict[str, Any]] = []
    aromatic_rows: List[Dict[str, Any]] = []
    ledger_rows: List[Dict[str, Any]] = []
    seen_for_ion: Set[str] = set()
    sources: Dict[str, Chem.Mol] = {}

    for operator in sorted(by_operator):
        available = by_operator[operator]
        ok, why = prefix_check(seed, ion.inventory_id, operator, available)
        if not ok:
            problems.append("ion %d, %s: %s" % (ion.inventory_id, operator,
                                                why))
        chosen = order_for(seed, ion.inventory_id, operator, available,
                           allowance)
        seen_for_operator: Set[str] = set()
        certified_count = 0
        distinct: Set[str] = set()

        for position, definition in enumerate(chosen, start=1):
            product = realise(labelled, definition)
            definition_id = "%d-%s-%d" % (ion.inventory_id, operator, position)
            examination_id = first_id + len(rows)
            row = {
                "examination_id": examination_id,
                "inventory_id": ion.inventory_id, "family_id": ion.family_id,
                "core_pattern": core_pattern, "grammar_id": grammar_id,
                "operator": operator, "mode": definition.mode,
                "order_index": position, "site": definition.site,
                "definition_id": definition_id,
                "parameters": definition.signature_text,
                "heavy_atom_change": definition.heavy_atom_change,
                "parent_smiles": parent_smiles,
                "parent_selfies": parent_selfies,
                "description": definition.description,
                "hydrogens_written_explicitly":
                    " ".join(str(i) for i in
                             product.hydrogens_written_explicitly),
            }

            # The ledger's working, written out. The conditions were decided
            # on these numbers before the product existed, so a reader who
            # cannot see the subtraction has to take the claim on trust.
            for result in definition.results:
                spent_after = (result.environment[3] - result.hydrogens)
                if result.index >= 0:
                    was = ion.mol.GetAtomWithIdx(result.index)
                    before_orders = " ".join(
                        sorted(str(b.GetBondType()) for b in was.GetBonds()))
                    hydrogens_before = int(was.GetTotalNumHs())
                    total_before = int(was.GetTotalValence())
                    spent_before = total_before - hydrogens_before
                    # Three roles, not two. An atom whose element the edit
                    # re-specifies is kept in place but is not the atom it
                    # was, so its total valence is not conserved and it
                    # answers to the creatable set like anything constructed.
                    role = ("element re-specified" if result.constructed
                            else "retained and touched")
                    atom = str(result.index)
                else:
                    before_orders, hydrogens_before = "", ""
                    total_before, spent_before = "", ""
                    role, atom = "constructed by the edit", "new"
                ledger_rows.append({
                    "examination_id": examination_id,
                    "inventory_id": ion.inventory_id,
                    "family_id": ion.family_id, "operator": operator,
                    "mode": definition.mode, "site": definition.site,
                    "atom": atom, "role": role,
                    "element": PT.GetElementSymbol(result.environment[0]),
                    "formal_charge": result.environment[1],
                    "aromatic": int(result.environment[2]),
                    "bond_orders_before": before_orders,
                    "hydrogens_before": hydrogens_before,
                    "valence_spent_before": spent_before,
                    "bond_orders_after": " ".join(result.environment[5]),
                    "hydrogens_after": result.hydrogens,
                    "valence_spent_after": spent_after,
                    "valence_spent_change": (
                        spent_after - spent_before
                        if spent_before != "" else ""),
                    "total_valence_before": total_before,
                    "total_valence_after": result.environment[3],
                    "answers_to": ("the creatable set" if result.constructed
                                   else "the observed set"),
                    "environment_after": "%s charge %+d%s valence %d %dH "
                                         "bonds %s"
                                         % (PT.GetElementSymbol(
                                                result.environment[0]),
                                            result.environment[1],
                                            " aromatic"
                                            if result.environment[2] else "",
                                            result.environment[3],
                                            result.hydrogens,
                                            " ".join(result.environment[5]))})

            if not product.built:
                row.update({
                    "certified": 0, "failing_requirement": product.failure,
                    "failing_reason": product.reason, "product_smiles": "",
                    "product_selfies": "", "mapped_product_smiles": "",
                    "identical_to_parent": "", "repeat_for_this_operator": "",
                    "repeat_for_this_ion": "", "similarity_to_parent": "",
                    "judge_accepted": "", "judge_reason": "",
                    "judge_agrees_with_certification": ""})
                rows.append(row)
                continue

            found = certify(product.mol, labelled, required_charge, expected,
                            len(ion.core))
            identical = product.smiles == parent_smiles
            verdict = judge(product.smiles, parent_smiles, query,
                            required_charge, scoring)
            # Certification reads the label applied before the edit; the judge
            # searches the product for the core. They answer the same question
            # by different routes, so a disagreement is worth recording.
            agrees = int(bool(found.certified)
                         == bool(verdict.core_retained
                                 and verdict.charge_on_core
                                 and verdict.required_charge
                                 and verdict.connected
                                 and verdict.closed_shell))
            row.update({
                "certified": int(found.certified),
                "failing_requirement": found.failed,
                "failing_reason": found.reason,
                "product_smiles": product.smiles,
                "product_selfies": selfies_of(product.smiles),
                "mapped_product_smiles": product.mapped_smiles,
                "identical_to_parent": int(identical),
                "repeat_for_this_operator":
                    int(product.smiles in seen_for_operator),
                "repeat_for_this_ion": int(product.smiles in seen_for_ion),
                "similarity_to_parent": similarity(parent_mol, product.mol),
                "judge_accepted": int(verdict.accepted),
                "judge_reason": verdict.reason,
                "judge_agrees_with_certification": agrees})
            rows.append(row)

            if found.certified:
                for attachment in definition.attachments:
                    if attachment.structure not in sources:
                        piece = Chem.MolFromSmiles(attachment.structure)
                        if piece is not None:
                            sources[attachment.structure] = piece
                tag = {"inventory_id": ion.inventory_id,
                       "family_id": ion.family_id,
                       "operator": operator, "mode": definition.mode,
                       "examination_id": examination_id}
                found_rings, found_aromatic = perceive(
                    product.mol, labelled, sources, tag)
                ring_rows.extend(found_rings)
                aromatic_rows.extend(found_aromatic)

            if found.certified and not identical:
                certified_count += 1
                distinct.add(product.smiles)
            seen_for_operator.add(product.smiles)
            seen_for_ion.add(product.smiles)

        summary[operator] = {
            "inventory_id": ion.inventory_id, "family_id": ion.family_id,
            "operator": operator, "definitions_available": len(available),
            "attempts_offered": len(chosen),
            "allowance": allowance,
            "examinations": len(chosen),
            "certified": certified_count,
            "distinct_products": len(distinct),
            "reached": int(certified_count > 0)}

    return rows, summary, problems, ring_rows, aromatic_rows, ledger_rows
