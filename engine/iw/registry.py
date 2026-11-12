
# -*- coding: utf-8 -*-
"""Stage 01b: the observed registries.

Both registries are read from the whole curated inventory, including ions that
end up without a core. An empty registry is an error rather than an open gate.

The creatable set is a subset of the observed set, not a separate valence
table: a state no curated ion displays is not reachable.

Rings are not edited, so no ring condition applies to an edit. The ring table
describes the inventory and is used by the perception check, which asserts that
no ring in a product differs from a ring in its parent.
"""
from __future__ import annotations

import collections
from typing import Any, Dict, List, Sequence, Tuple

from rdkit import Chem

from .common import (AtomEnv, BondEnv, RingEnv, RunError, atom_environment,
                     bond_environment, bond_row, env_to_row,
                     parse, ring_environments)

PT = Chem.GetPeriodicTable()

HALOGENS = frozenset((9, 17, 35, 53, 85))

ATOM_ENV_COLUMNS = (
    "environment_id", "element", "atomic_number", "formal_charge", "aromatic",
    "total_valence", "attached_hydrogens", "bond_orders", "degree",
    "occurrences", "ions", "creatable", "not_creatable_reason",
)

BOND_ENV_COLUMNS = (
    "bond_environment_id", "first_element", "first_formal_charge",
    "first_aromatic", "second_element", "second_formal_charge",
    "second_aromatic", "bond_order", "occurrences", "ions",
)

RING_ENV_COLUMNS = ("bond_order", "ring_size", "occurrences", "ions")


def lowest_standard_valence(atomic_number: int) -> int:
    values = [v for v in PT.GetValenceList(atomic_number) if v > 0]
    return min(values) if values else PT.GetDefaultValence(atomic_number)


#: The conditions a placed atom's environment must clear, in the order they are
#: applied. The definition space surviving each one is reported, so the
#: order is part of the method and the codes are what the shrinkage table keys
#: on. `observed` is not in this list: it is the atom-environment condition
#: itself, which needs the registry and is applied by `placement`.
RESTRICTIONS: Tuple[str, ...] = ("neutral_and_non_aromatic",
                                 "terminal_halogen",
                                 "above_lowest_valence_with_hydrogen")

REASONS: Dict[str, str] = {
    "charged": "carries a formal charge; only the core carries charge",
    "aromatic": "aromatic; an edit places no aromatic atom",
    "terminal_halogen": "neutral halogen away from a terminal single bond",
    "above_lowest_valence_with_hydrogen":
        "above the lowest standard valence while carrying hydrogen; "
        "sanitization can rewrite such a centre with its charge separated",
    "not_observed": "no curated ion displays this atom environment",
    "bond_not_observed": "no curated ion displays a bond between these "
                         "two kinds of atom at this order",
}


def restriction_failed(env: AtomEnv) -> Tuple[str, str]:
    """The first restriction this environment fails, and why.

    Returned as a code so the shrinkage table of C4 can say which condition did
    the work, instead of reporting one refusal count for all of them.
    """
    z, charge, aromatic, valence, hydrogens, orders = env

    if charge != 0:
        return "neutral_and_non_aromatic", REASONS["charged"]
    if aromatic:
        return "neutral_and_non_aromatic", REASONS["aromatic"]
    if z in HALOGENS and (len(orders) != 1 or orders[0] != "SINGLE"):
        return "terminal_halogen", REASONS["terminal_halogen"]
    if valence > lowest_standard_valence(z) and hydrogens > 0:
        return ("above_lowest_valence_with_hydrogen",
                REASONS["above_lowest_valence_with_hydrogen"])
    return "", ""


def creatable(env: AtomEnv) -> Tuple[bool, str]:
    """Whether an edit may place an atom into this observed environment.

    Only neutral, non-aromatic atoms are placed, and two further
    restrictions apply, each with its cause. One implementation, read two
    ways: this answers yes or no, and `restriction_failed` says which
    condition said no.
    """
    code, reason = restriction_failed(env)
    return (not code), reason


def build(inventory: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Both registries, from every curated ion."""
    if not inventory:
        raise RunError("the registries are learned from the inventory and the "
                       "inventory is empty; an unlearned registry is an error")

    atom_counts: Dict[AtomEnv, int] = collections.Counter()
    atom_ions: Dict[AtomEnv, set] = collections.defaultdict(set)
    ring_counts: Dict[RingEnv, int] = collections.Counter()
    ring_ions: Dict[RingEnv, set] = collections.defaultdict(set)
    bond_counts: Dict[BondEnv, int] = collections.Counter()
    bond_ions: Dict[BondEnv, set] = collections.defaultdict(set)

    for record in inventory:
        iid = record["inventory_id"]
        mol = parse(record["canonical_smiles"])
        for atom in mol.GetAtoms():
            env = atom_environment(atom)
            atom_counts[env] += 1
            atom_ions[env].add(iid)
        for bond in mol.GetBonds():
            benv = bond_environment(bond)
            bond_counts[benv] += 1
            bond_ions[benv].add(iid)
        for renv in ring_environments(mol):
            ring_counts[renv] += 1
            ring_ions[renv].add(iid)

    if not atom_counts:
        raise RunError("no atom environments were learned")

    atom_rows: List[Dict[str, Any]] = []
    for index, env in enumerate(sorted(atom_counts), start=1):
        ok, reason = creatable(env)
        row = env_to_row(env)
        row.update({
            "environment_id": index,
            "occurrences": atom_counts[env],
            "ions": len(atom_ions[env]),
            "creatable": int(ok),
            "not_creatable_reason": reason,
        })
        atom_rows.append(row)

    bond_rows: List[Dict[str, Any]] = []
    for index, env in enumerate(sorted(bond_counts), start=1):
        row = bond_row(env)
        row.update({
            "bond_environment_id": index,
            "occurrences": bond_counts[env],
            "ions": len(bond_ions[env]),
        })
        bond_rows.append(row)

    ring_rows = [{
        "bond_order": order,
        "ring_size": size,
        "occurrences": ring_counts[(order, size)],
        "ions": len(ring_ions[(order, size)]),
    } for order, size in sorted(ring_counts)]

    creatable_set = {env for env in atom_counts if creatable(env)[0]}

    return {
        "atom_rows": atom_rows,
        "bond_rows": bond_rows,
        "ring_rows": ring_rows,
        "observed": set(atom_counts),
        "observed_bonds": set(bond_counts),
        "creatable": creatable_set,
        "observed_rings": set(ring_counts),
        "summary": {
            "ions_read": len(inventory),
            "atom_environments_observed": len(atom_counts),
            "atom_environments_creatable": len(creatable_set),
            "bond_environments_observed": len(bond_counts),
            "ring_environments_observed": len(ring_counts),
        },
    }

