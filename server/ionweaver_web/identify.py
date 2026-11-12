"""Recognise a pasted ion by the 44 published cores.

The ion is put through the same per-structure checks curation applies to the
source inventory (one connected, closed-shell, metal-free structure with a net
charge of +1 or −1 carried by one atom).  Every published core of its ion type is
then matched with the engine's own ``cores.rank_matches``, which accepts a
match only when all formal charge lies on the matched core and everything
outside it is neutral.  Matches are ranked by core size, then by support: a
larger core fixes more of the ion and is the more specific statement of family
membership.  An ion that is itself in the reference inventory is offered under
the core the published run assigned to it.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from rdkit import Chem

from iw import cores as C
from iw import fragments
from iw.common import CHEMISTRY_ERRORS, canonical, fragment_count, has_metal, net_charge, parse, radical_electrons

from .library import Published

__all__ = ["check_ion", "identify", "make_record", "REASONS"]

REASONS = {
    "unreadable_structure": "The SMILES string could not be read as a molecule.",
    "disconnected": "The structure has more than one disconnected part. Give one ion, without its counter-ion.",
    "radical": "The structure carries an unpaired electron.",
    "metal_atom": "The structure contains a metal atom; the method edits organic ions only.",
    "not_monovalent": "The net charge is not +1 or −1; the method covers monovalent ions.",
    "charge_not_on_a_single_atom": "The charge is spread over more than one atom (for example a nitro group written with two charged atoms); the method holds one charged atom.",
}


def check_ion(smiles: str) -> Dict[str, Any]:
    """The per-structure checks of curation, on one pasted ion."""
    text = (smiles or "").strip()
    if not text:
        return {"ok": False, "status": "unreadable_structure", "reason": "Paste a SMILES string.", "input": text}
    try:
        mol = parse(text)
    except CHEMISTRY_ERRORS:
        return {"ok": False, "status": "unreadable_structure", "reason": REASONS["unreadable_structure"], "input": text}
    if fragment_count(mol) != 1:
        return {"ok": False, "status": "disconnected", "reason": REASONS["disconnected"], "input": text}
    if radical_electrons(mol) != 0:
        return {"ok": False, "status": "radical", "reason": REASONS["radical"], "input": text}
    if has_metal(mol):
        return {"ok": False, "status": "metal_atom", "reason": REASONS["metal_atom"], "input": text}
    charge = net_charge(mol)
    if charge not in (1, -1):
        return {"ok": False, "status": "not_monovalent", "reason": REASONS["not_monovalent"], "input": text,
                "charge": charge}
    if sum(1 for a in mol.GetAtoms() if int(a.GetFormalCharge()) != 0) != 1:
        return {"ok": False, "status": "charge_not_on_a_single_atom",
                "reason": REASONS["charge_not_on_a_single_atom"], "input": text, "charge": charge}
    canon = canonical(mol)
    return {"ok": True, "status": "retained", "reason": "", "input": text, "canonical_smiles": canon,
            "charge": charge, "ion_type": "cation" if charge > 0 else "anion", "mol": parse(canon)}


def identify(smiles: str, P: Published) -> Dict[str, Any]:
    cur = check_ion(smiles)
    if not cur["ok"]:
        return {k: v for k, v in cur.items() if k != "mol"}
    mol = cur["mol"]
    matches: List[Dict[str, Any]] = []
    for cm in P.core_models:
        if cm.ion_type != cur["ion_type"] or cm.query is None:
            continue
        best, _trace = C.rank_matches(mol, cm.query, cur["charge"])
        if best is None:
            continue
        matches.append({**cm.as_dict(), "core_atoms": [int(i) for i in best]})
    out = {k: v for k, v in cur.items() if k != "mol"}
    if not matches:
        out.update({"ok": False, "status": "no_core",
                    "reason": "No published family core matches this ion with all of its charge on the core. "
                              "Upload a sheet with at least two related ions of the same family so IonWeaver "
                              "can discover the core they share."})
        return out
    matches.sort(key=lambda m: (-m["atoms"], -m["support"], m["family_label"], m["grammar_id"]))
    known = P.record_by_smiles(cur["canonical_smiles"])
    if known is not None:
        i = next((k for k, m in enumerate(matches)
                  if m["family_id"] == known["family_id"] and m["grammar_id"] == known["grammar_id"]), None)
        if i is not None:
            matches.insert(0, matches.pop(i))
    out.update({"best": matches[0], "alternatives": matches[1:],
                "published_id": known["inventory_id"] if known is not None else None})
    return out


def make_record(canonical_smiles: str, family_id: str, family_label: str, grammar_id: str,
                pattern: str, core_atoms: List[int], inventory_id: int = -1) -> Dict[str, Any]:
    """A starting-ion record of the shape ``fragments.starting_ions`` produces."""
    try:
        mol = parse(canonical_smiles)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("That SMILES could not be read.") from exc
    core = sorted({int(i) for i in core_atoms})
    if not core or core[0] < 0 or core[-1] >= mol.GetNumAtoms():
        raise ValueError("The core atoms do not index the ion's atoms.")
    if not str(family_id).startswith(("cation", "anion")):
        raise ValueError("The family must be a cation or an anion family.")
    if Chem.MolFromSmarts(pattern or "") is None:
        raise ValueError("The core pattern is not valid SMARTS.")
    fixed, absorbed = fragments.fixed_atoms(mol, core)
    pieces = fragments.components_of(mol, core)
    return {
        "inventory_id": int(inventory_id), "family_id": family_id, "family_label": family_label,
        "canonical_smiles": canonical_smiles, "grammar_id": grammar_id, "core_pattern": pattern,
        "core_atoms": " ".join(str(i) for i in core), "core_atom_count": len(core),
        "ring_shared_atoms": " ".join(str(i) for i in sorted(absorbed)), "ring_shared_count": len(absorbed),
        "heavy_atoms": int(mol.GetNumAtoms()), "body_heavy_atoms": int(mol.GetNumAtoms() - len(fixed)),
        "body_components": len(pieces),
    }
