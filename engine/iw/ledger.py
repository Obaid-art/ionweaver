
# -*- coding: utf-8 -*-
"""The valence ledger, and the atom environment an edit will produce.

Every condition in the method is decided here, from the parent graph
and the change an edit specifies, with no product built. An edit that cannot be
given a legal valence state is never created, and one whose touched atoms would
land outside the sets the inventory displays is never created either.

The ledger computes the hydrogen count and the builder writes it. Nothing
relies on a sanitization pass recomputing implicit hydrogens, because for a
charged or aromatic atom RDKit does not: it keeps the explicit count it was
given and the product goes over-valent instead. Predicting the count is what
makes the condition decidable in advance, and the builder then applies it, so
the two can never disagree silently. After a product is built its touched atoms
are compared against these predictions, and a mismatch stops the run.

Total valence is conserved at every retained atom. The freed or consumed valence
is taken up by hydrogen, so an edit changes how an atom's valence is spent and
never how much of it there is. That is why a single arithmetic rule covers every
operator, and why an edit whose arithmetic needs a negative hydrogen count is
refused before anything is assembled.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

from .common import AtomEnv, RunError

PT = Chem.GetPeriodicTable()

#: How much valence each bond order spends. AROMATIC is absent on purpose: no
#: operator alters a ring bond, so no aromatic bond is ever added or removed,
#: and an aromatic order reaching the arithmetic is a defect that stops the run.
ORDER_VALUE: Dict[str, int] = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}

BOND_TYPE: Dict[str, Chem.BondType] = {
    "SINGLE": Chem.BondType.SINGLE,
    "DOUBLE": Chem.BondType.DOUBLE,
    "TRIPLE": Chem.BondType.TRIPLE,
}


def order_value(order: str) -> int:
    if order not in ORDER_VALUE:
        raise RunError(
            "the valence ledger was handed a bond of order %r. Only single, "
            "double and triple bonds are ever added or removed, because no "
            "operator alters a ring bond." % order)
    return ORDER_VALUE[order]


def lowest_valence_at_or_above(atomic_number: int, spent: int) -> Optional[int]:
    """The valence RDKit will settle an atom into, given what its bonds spend.

    The periodic table's valence list, in order, is what an implicit hydrogen
    count is computed from: the first allowed valence that covers the bonds.
    A list entry of -1 means the element accepts any valence, and then nothing
    is added.
    """
    for value in PT.GetValenceList(atomic_number):
        if value < 0:
            return spent
        if value >= spent:
            return value
    return None


@dataclass
class Touch:
    """What an edit does to one atom it keeps.

    `removed` and `added` are the orders of the bonds at this atom that the edit
    takes away and puts back. `new_element` is set by Substitute alone, which
    keeps every bond and changes only what the atom is.
    """
    index: int
    removed: Tuple[str, ...] = ()
    added: Tuple[str, ...] = ()
    new_element: Optional[int] = None

    @property
    def spent_delta(self) -> int:
        return (sum(order_value(o) for o in self.added)
                - sum(order_value(o) for o in self.removed))


@dataclass
class Placed:
    """An atom the edit creates, described before it exists."""
    atomic_number: int
    bonds: Tuple[str, ...]
    label: str = ""


@dataclass
class Resulting:
    """One atom's predicted state, and the hydrogen count the builder applies.

    `constructed` says which set the environment is tested against. An atom the
    edit builds, and an atom whose element it re-specifies, are both constructed
    and answer to the creatable set. An atom it keeps as it is, and merely
    touches, answers to the observed set, because it was already in whatever
    state it was in.
    """
    index: int
    environment: AtomEnv
    hydrogens: int
    constructed: bool
    legal: bool = True
    reason: str = ""


def _orders_after(atom: Chem.Atom, touch: Touch) -> Optional[Tuple[str, ...]]:
    """The atom's bond orders once the edit has run."""
    orders: List[str] = sorted(str(b.GetBondType()) for b in atom.GetBonds())
    for order in touch.removed:
        if order not in orders:
            return None
        orders.remove(order)
    orders.extend(touch.added)
    return tuple(sorted(orders))


def retained(mol: Chem.Mol, touch: Touch) -> Resulting:
    """The environment a retained, touched atom will be left in.

    Element, charge and aromatic status carry over, except that Substitute sets
    the element. Valence is conserved and hydrogen absorbs the difference, so
    the whole prediction is one subtraction. A negative hydrogen count means the
    edit has no legal valence state at this atom, and it is never created.
    """
    atom = mol.GetAtomWithIdx(int(touch.index))
    orders = _orders_after(atom, touch)
    if orders is None:
        raise RunError(
            "the ledger was told to remove a bond of order %s from atom %d of "
            "%s, which does not carry one" % (touch.removed, touch.index,
                                              Chem.MolToSmiles(mol)))

    charge = int(atom.GetFormalCharge())
    aromatic = bool(atom.GetIsAromatic())
    hydrogens_before = int(atom.GetTotalNumHs())
    spent_before = int(atom.GetTotalValence()) - hydrogens_before

    if touch.new_element is not None:
        # Substitute keeps every bond, so what the bonds spend is unchanged and
        # the hydrogen count is the one the new element's valence list gives.
        element = int(touch.new_element)
        valence = lowest_valence_at_or_above(element, spent_before)
        if valence is None:
            return Resulting(touch.index, (element, charge, aromatic, 0, 0,
                                           orders), 0, True, False,
                             "%s has no standard valence that covers the %d "
                             "its bonds spend"
                             % (PT.GetElementSymbol(element), spent_before))
        hydrogens = valence - spent_before
        return Resulting(touch.index,
                         (element, charge, aromatic, valence, hydrogens,
                          orders), hydrogens, True)

    element = int(atom.GetAtomicNum())
    hydrogens = hydrogens_before - touch.spent_delta
    if hydrogens < 0:
        return Resulting(
            touch.index,
            (element, charge, aromatic, int(atom.GetTotalValence()),
             hydrogens, orders), hydrogens, False, False,
            "the edit spends %d more valence at this atom than its %d hydrogen"
            "(s) can give up" % (touch.spent_delta, hydrogens_before))
    valence = spent_before + touch.spent_delta + hydrogens
    return Resulting(touch.index,
                     (element, charge, aromatic, valence, hydrogens, orders),
                     hydrogens, False)


def placed(atom: Placed, index: int = -1) -> Resulting:
    """The environment an atom the edit creates will be left in.

    A placed atom is neutral and non-aromatic, per B4, so its hydrogen count
    follows from its element and from what its bonds spend.
    """
    spent = sum(order_value(o) for o in atom.bonds)
    valence = lowest_valence_at_or_above(atom.atomic_number, spent)
    orders = tuple(sorted(atom.bonds))
    if valence is None:
        return Resulting(index, (atom.atomic_number, 0, False, 0, 0, orders),
                         0, True, False,
                         "%s has no standard valence that covers the %d its "
                         "bonds would spend"
                         % (PT.GetElementSymbol(atom.atomic_number), spent))
    hydrogens = valence - spent
    return Resulting(index,
                     (atom.atomic_number, 0, False, valence, hydrogens,
                      orders), hydrogens, True)
