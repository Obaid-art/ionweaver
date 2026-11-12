
# -*- coding: utf-8 -*-
"""The ring and aromaticity invariant, re-perceived on every product.

Two properties follow from the operators themselves rather than from any
condition, and both are asserted on finished products:

The ring invariant. Every ring in a product is a ring of the parent or a
ring of the fragment the edit introduced, unchanged in size and in bond orders.
No ring is built and none is altered. Because every library component is cut
from a curated ion, every ring in every product traces back to a curated ion.

The aromaticity invariant. No atom's aromatic flag differs from the flag it
carried in the parent, or in the fragment introduced.

Nothing is restricted pre-emptively. One path cannot be ruled out by
argument alone: a change to an exocyclic bond on a ring atom. No restriction is
written for it. If this check fires, it names the operator, the ion and the atom
that needs one, and the restriction is added then, with its cost in definitions
reported. No reach is surrendered to a hypothesis.

A non-empty mismatch here is a failed run and not a finding.
"""
from __future__ import annotations

from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

from .build import ADDED, CORE, PARENT, SOURCE

RING_COLUMNS = (
    "inventory_id", "family_id", "operator", "mode", "examination_id",
    "ring_index", "ring_size", "bond_orders", "source", "source_ring_size",
    "source_bond_orders", "agrees", "note")

AROMATIC_COLUMNS = (
    "inventory_id", "family_id", "operator", "mode", "examination_id",
    "atom", "element", "source", "aromatic_in_product",
    "aromatic_in_source", "agrees")


def _ring_signature(mol: Chem.Mol, ring: Sequence[int]) -> Tuple[int, str]:
    """A ring's size and the multiset of its bond orders."""
    members = list(ring)
    orders: List[str] = []
    for position, atom in enumerate(members):
        other = members[(position + 1) % len(members)]
        bond = mol.GetBondBetweenAtoms(int(atom), int(other))
        if bond is None:
            return len(members), "the ring is not closed in this molecule"
        orders.append(str(bond.GetBondType()))
    return len(members), " ".join(sorted(orders))


def _rings_by_atoms(mol: Chem.Mol) -> Dict[FrozenSet[int], Sequence[int]]:
    return {frozenset(int(i) for i in ring): ring
            for ring in mol.GetRingInfo().AtomRings()}


def check(product: Chem.Mol, parent: Chem.Mol,
          sources: Dict[str, Chem.Mol],
          tag: Dict[str, Any]) -> Tuple[List[Dict[str, Any]],
                                        List[Dict[str, Any]]]:
    """One product's rings and aromatic flags, each beside its source's.

    `sources` maps the structure string an edit attached to that structure as a
    molecule, so an atom the edit introduced can be compared with the atom it
    was copied from.
    """
    ring_rows: List[Dict[str, Any]] = []
    aromatic_rows: List[Dict[str, Any]] = []

    parent_rings = {frozenset(int(i) for i in ring): ring
                    for ring in parent.GetRingInfo().AtomRings()}
    source_rings = {key: _rings_by_atoms(mol) for key, mol in sources.items()}

    for number, ring in enumerate(product.GetRingInfo().AtomRings(), start=1):
        size, orders = _ring_signature(product, ring)
        labels: List[int] = []
        added: List[Tuple[str, int]] = []
        for index in ring:
            atom = product.GetAtomWithIdx(int(index))
            if atom.HasProp(PARENT):
                labels.append(int(atom.GetIntProp(PARENT)))
            elif atom.HasProp(ADDED):
                added.append((atom.GetProp(SOURCE),
                              int(atom.GetIntProp(ADDED))))

        row = dict(tag, ring_index=number, ring_size=size,
                   bond_orders=orders, source="", source_ring_size="",
                   source_bond_orders="", agrees=0, note="")

        if labels and added:
            row["source"] = "both the parent and an introduced fragment"
            row["note"] = ("this ring spans atoms of the parent and atoms the "
                           "edit introduced, so the edit created a ring")
            ring_rows.append(row)
            continue

        if labels:
            key = frozenset(labels)
            found = parent_rings.get(key)
            row["source"] = "parent"
            if found is None:
                row["note"] = ("the same atoms of the parent do not form a "
                               "ring there")
                ring_rows.append(row)
                continue
            was_size, was_orders = _ring_signature(parent, found)
            row["source_ring_size"] = was_size
            row["source_bond_orders"] = was_orders
            row["agrees"] = int(was_size == size and was_orders == orders)
            if not row["agrees"]:
                row["note"] = "the ring's size or its bond orders changed"
            ring_rows.append(row)
            continue

        structures = {s for s, _ in added}
        if len(structures) != 1:
            row["source"] = "more than one introduced fragment"
            row["note"] = ("this ring spans atoms from %d introduced "
                           "structures" % len(structures))
            ring_rows.append(row)
            continue
        structure = next(iter(structures))
        key = frozenset(i for _, i in added)
        found = source_rings.get(structure, {}).get(key)
        row["source"] = "the fragment %s" % structure
        if found is None:
            row["note"] = ("the same atoms of %s do not form a ring there"
                           % structure)
            ring_rows.append(row)
            continue
        was_size, was_orders = _ring_signature(sources[structure], found)
        row["source_ring_size"] = was_size
        row["source_bond_orders"] = was_orders
        row["agrees"] = int(was_size == size and was_orders == orders)
        if not row["agrees"]:
            row["note"] = "the ring's size or its bond orders changed"
        ring_rows.append(row)

    for atom in product.GetAtoms():
        aromatic = int(bool(atom.GetIsAromatic()))
        if atom.HasProp(PARENT):
            label = int(atom.GetIntProp(PARENT))
            was = int(bool(parent.GetAtomWithIdx(label).GetIsAromatic()))
            where = "parent atom %d" % label
        elif atom.HasProp(ADDED):
            structure = atom.GetProp(SOURCE)
            index = int(atom.GetIntProp(ADDED))
            piece = sources.get(structure)
            if piece is None:
                was, where = -1, "the fragment %s, which was not supplied" \
                    % structure
            else:
                was = int(bool(piece.GetAtomWithIdx(index).GetIsAromatic()))
                where = "%s atom %d" % (structure, index)
        else:
            was, where = -1, "no source"
        if was == aromatic:
            continue        # agreement is the expected case and is not written
        aromatic_rows.append(dict(
            tag, atom=int(atom.GetIdx()), element=atom.GetSymbol(),
            source=where, aromatic_in_product=aromatic,
            aromatic_in_source=was, agrees=0))

    return ring_rows, aromatic_rows


def mismatches(ring_rows: Sequence[Dict[str, Any]],
               aromatic_rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Every row that did not agree. This must be empty."""
    out = [row for row in ring_rows if not row["agrees"]]
    out.extend(row for row in aromatic_rows if not row["agrees"])
    return out
