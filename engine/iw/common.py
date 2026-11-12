
# -*- coding: utf-8 -*-
"""Shared primitives, imported by every stage.

Every reported quantity is computed from the records rather than carried as a
literal.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# Spec Part A: the valence, aromaticity and ring-perception models are RDKit's.
# A different version is a different experiment, so the version is recorded in
# every settings file and in the manifest.
RDKIT_VERSION = Chem.rdBase.rdkitVersion

METALS = frozenset(
    list(range(3, 5)) + list(range(11, 14)) + list(range(19, 32))
    + list(range(37, 51)) + list(range(55, 85)) + list(range(87, 104))
)

# Chemistry errors that may be recorded. Anything else stops the run.
CHEMISTRY_ERRORS = (ValueError, Chem.rdchem.MolSanitizeException)


class RunError(RuntimeError):
    """A failure that must stop the run rather than become a refusal record."""


# ---------------------------------------------------------------- molecules

def parse(smiles: str) -> Chem.Mol:
    """Read and sanitize. Raises a chemistry error if the text is not a molecule."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("unreadable_smiles")
    return mol


def normalise(mol: Chem.Mol) -> Chem.Mol:
    """Strip stereochemistry and isotope labels, in place, before canonicalising.

    Neither is chemistry this method has any semantics for. No operator sets a
    stereocentre or an isotope, the atom-environment condition keys on atomic
    number, and a family of ionic liquids is not distinguished by either. Left
    in, a deuterated analogue and its protiated form become two starting ions
    with the same molecular graph, and two enantiomers become two ions whose
    products are indistinguishable to every condition the engine applies.
    """
    Chem.RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    # Clearing the isotope leaves a former [2H] as an explicit hydrogen node,
    # so it has to be folded back into its neighbour's implicit count. Without
    # this a deuterated analogue stays a separate graph under a different
    # spelling, which is the opposite of what stripping the label is for.
    return Chem.RemoveHs(mol)


def canonical(mol: Chem.Mol) -> str:
    """The canonical form the whole study keys on, after normalisation."""
    return Chem.MolToSmiles(normalise(Chem.Mol(mol)))


def net_charge(mol: Chem.Mol) -> int:
    return sum(int(a.GetFormalCharge()) for a in mol.GetAtoms())


def radical_electrons(mol: Chem.Mol) -> int:
    return sum(int(a.GetNumRadicalElectrons()) for a in mol.GetAtoms())


def fragment_count(mol: Chem.Mol) -> int:
    return len(Chem.GetMolFrags(mol))


def has_metal(mol: Chem.Mol) -> bool:
    return any(int(a.GetAtomicNum()) in METALS for a in mol.GetAtoms())


def heavy_atoms(mol: Chem.Mol) -> int:
    return mol.GetNumHeavyAtoms()


# ------------------------------------------------- atom and ring environments

#: An atom environment, as spec C1 defines it. Hashable, and the only key the
#: creatable set is ever derived from.
AtomEnv = Tuple[int, int, bool, int, int, Tuple[str, ...]]


def atom_environment(atom: Chem.Atom) -> AtomEnv:
    """Element, formal charge, aromatic status, total valence, attached
    hydrogens, and the ordered multiset of bond orders this atom carries."""
    orders = tuple(sorted(str(b.GetBondType()) for b in atom.GetBonds()))
    return (
        int(atom.GetAtomicNum()),
        int(atom.GetFormalCharge()),
        bool(atom.GetIsAromatic()),
        int(atom.GetTotalValence()),
        int(atom.GetTotalNumHs()),
        orders,
    )


#: A bond environment, radius one: what each end is, and the order between
#: them. An end is the element, the formal charge and the aromatic flag, and
#: nothing else. Hydrogens and the rest of an end's bonds are already governed
#: by the atom-environment condition, so repeating them here would only make
#: the same test twice and narrow the reachable set for no further reason.
BondEnd = Tuple[int, int, bool]
BondEnv = Tuple[BondEnd, BondEnd, str]


def bond_end(atom: Chem.Atom) -> BondEnd:
    return (int(atom.GetAtomicNum()), int(atom.GetFormalCharge()),
            bool(atom.GetIsAromatic()))


def bond_key(first: BondEnd, second: BondEnd, order: str) -> BondEnv:
    """The pair, written the same way whichever end is given first."""
    low, high = sorted((first, second))
    return (low, high, order)


def bond_environment(bond: Chem.Bond) -> BondEnv:
    return bond_key(bond_end(bond.GetBeginAtom()),
                    bond_end(bond.GetEndAtom()),
                    str(bond.GetBondType()))


def bond_row(env: BondEnv) -> Dict[str, Any]:
    table = Chem.GetPeriodicTable()
    first, second, order = env
    return {
        "first_element": table.GetElementSymbol(first[0]),
        "first_formal_charge": first[1],
        "first_aromatic": int(first[2]),
        "second_element": table.GetElementSymbol(second[0]),
        "second_formal_charge": second[1],
        "second_aromatic": int(second[2]),
        "bond_order": order,
    }


def env_to_row(env: AtomEnv) -> Dict[str, Any]:
    z, chg, arom, val, hs, orders = env
    return {
        "atomic_number": z,
        "element": Chem.GetPeriodicTable().GetElementSymbol(z),
        "formal_charge": chg,
        "aromatic": int(arom),
        "total_valence": val,
        "attached_hydrogens": hs,
        "bond_orders": " ".join(orders),
        "degree": len(orders),
    }


#: A ring environment: one bond order together with the size of the ring the
#: bond sits in. The pair is required, never the smallest ring of a
#: product's SSSR.
RingEnv = Tuple[str, int]


def ring_environments(mol: Chem.Mol) -> List[RingEnv]:
    """Every (bond order, ring size) pair the molecule displays."""
    out: List[RingEnv] = []
    info = mol.GetRingInfo()
    for ring in info.AtomRings():
        size = len(ring)
        members = set(ring)
        for bond in mol.GetBonds():
            if (bond.GetBeginAtomIdx() in members
                    and bond.GetEndAtomIdx() in members
                    and bond.IsInRingSize(size)):
                out.append((str(bond.GetBondType()), size))
    return out


# ----------------------------------------------------------------- the judge

@dataclass(frozen=True)
class Scoring:
    """The scoring settings, passed explicitly at every scoring site.

    The fields have no defaults, so a divergence between two sites is visible
    at the call rather than implicit in it.
    """
    require_closed_shell: bool
    require_connected: bool
    exclude_parent: bool
    allowance: Optional[int]


#: The only instance the study uses. Every stage imports this; nothing
#: constructs its own.
STUDY_SCORING = Scoring(require_closed_shell=True, require_connected=True,
                        exclude_parent=True, allowance=50)


@dataclass(frozen=True)
class Verdict:
    sanitises: bool
    connected: bool
    closed_shell: bool
    required_charge: bool
    core_retained: bool
    charge_on_core: bool
    identical_to_parent: bool
    accepted: bool
    reason: str = ""

    def as_row(self) -> Dict[str, Any]:
        return {
            "sanitises": int(self.sanitises),
            "connected": int(self.connected),
            "closed_shell": int(self.closed_shell),
            "required_charge": int(self.required_charge),
            "core_retained": int(self.core_retained),
            "charge_on_core": int(self.charge_on_core),
            "identical_to_parent": int(self.identical_to_parent),
            "accepted": int(self.accepted),
            "reason": self.reason,
        }


def judge(product_smiles: str,
          parent_canonical: str,
          core_pattern: Optional[Chem.Mol],
          required_charge: int,
          scoring: Scoring) -> Verdict:
    """The single scoring site. Every method is judged here and nowhere else.

    `core_pattern` is a SMARTS query. Passing None tests everything except the
    two core conditions, which then report True, so a caller that has no core
    cannot silently acquire a pass it did not earn -- it must not ask.
    """
    try:
        mol = parse(product_smiles)
    except CHEMISTRY_ERRORS:
        return Verdict(False, False, False, False, False, False, False, False,
                       "does_not_sanitise")

    connected = fragment_count(mol) == 1
    closed_shell = radical_electrons(mol) == 0
    charge_ok = net_charge(mol) == int(required_charge)
    identical = canonical(mol) == str(parent_canonical)

    core_ok, loc_ok = True, True
    if core_pattern is not None:
        match = _best_core_match(mol, core_pattern)
        core_ok = match is not None
        if core_ok:
            charged = {a.GetIdx() for a in mol.GetAtoms()
                       if int(a.GetFormalCharge()) != 0}
            loc_ok = bool(charged) and charged <= set(match)
        else:
            loc_ok = False

    failures: List[str] = []
    if scoring.require_connected and not connected:
        failures.append("disconnected")
    if scoring.require_closed_shell and not closed_shell:
        failures.append("radical")
    if not charge_ok:
        failures.append("charge")
    if not core_ok:
        failures.append("core_absent")
    elif not loc_ok:
        failures.append("charge_off_core")
    if scoring.exclude_parent and identical:
        failures.append("identical_to_parent")

    return Verdict(True, connected, closed_shell, charge_ok, core_ok, loc_ok,
                   identical, not failures, ",".join(failures))


def _best_core_match(mol: Chem.Mol, pattern: Chem.Mol):
    """Locate the core in a product by substructure match.

    This is the independent route, used by the judge. Certification does not
    come this way: requirements 3 and 4 read the label applied before the edit,
    so no fresh search is run on a product. Prefer a match carrying formal
    charge, so a pattern matching in more than one place cannot be scored
    against an uncharged copy of itself.
    """
    hits = mol.GetSubstructMatches(pattern, uniquify=False)
    if not hits:
        return None
    for h in hits:
        if any(int(mol.GetAtomWithIdx(i).GetFormalCharge()) for i in h):
            return h
    return hits[0]


# ------------------------------------------------------------------ plumbing

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: str, rows: Sequence[Dict[str, Any]],
              columns: Optional[Sequence[str]] = None) -> int:
    """Write rows and return the count. An empty file is a result, so a zero-row
    table is still written, with its header, whenever the columns are known."""
    import csv
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cols = list(columns) if columns else (list(rows[0].keys()) if rows else [])
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return len(rows)


def write_json(path: str, payload: Any) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
