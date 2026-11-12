
# -*- coding: utf-8 -*-
"""Stage 03 — core discovery.

Core discovery uses one bond comparator, one match ranking with
symmetry-equivalent matches collapsed, a component-wise charge rule outside the
core, two query forms compared by support, a floor of three atoms and two
bonds, two threshold ladders and a pairwise residual limit.

Each of the four tie-break decision points is recorded. The maximum common
substructure search runs with no wall-clock limit and asserts that it was never
cancelled, so a truncated search cannot pass unnoticed.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from .common import RunError

#: Exhaustive pairwise MCS is used when a family residual is at most this size.
#: It keeps discovery global and order-independent for small or heterogeneous
#: residuals; larger ones rely on the threshold and charge-environment passes
#: and reach this branch once they shrink.
PAIRWISE_RESIDUAL_LIMIT = 24

#: A candidate below this size is a functional group, not a substructure that
#: can anchor a family identity. Without this floor a bare heteroatom wins on
#: support alone: a one-atom boron beat a fluoroborate core, and a two-atom
#: acyclic fragment beat two triazolium ring cores.
MIN_CORE_ATOMS = 3
MIN_CORE_BONDS = 2

#: Support below this is a pairwise overlap rather than a family motif.
MIN_SUPPORT = 2


def candidate_thresholds(n: int) -> List[float]:
    """The support ladder for the whole-residual pass, with a floor that falls
    as the residual grows so two members can still carry a candidate."""
    values = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50,
              max(2.0 / max(n, 2), 0.10)]
    return sorted({round(min(1.0, max(0.01, x)), 4) for x in values},
                  reverse=True)


def charge_environment_thresholds(n: int) -> Tuple[float, ...]:
    """The ladder for the charge-environment pass. Shorter than the above
    because each group is already homogeneous in how its charge sits."""
    return (1.0, 0.8, max(2.0 / max(n, 1), 0.5))


def mcs_result(mols: Sequence[Chem.Mol], threshold: float) -> Optional[Chem.Mol]:
    """The shared substructure of `mols` at `threshold` support.

    No wall-clock limit is set. Under a limit, a truncated substructure
    could become a family's core with nothing recording it, and a slower
    machine would give a different answer.
    """
    if len(mols) < 2:
        return None
    params = rdFMCS.MCSParameters()
    params.Timeout = 0
    params.Threshold = float(threshold)
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrder
    params.AtomCompareParameters.MatchFormalCharge = True
    params.AtomCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    params.MaximizeBonds = True
    result = rdFMCS.FindMCS(list(mols), params)
    if result.canceled:
        raise RunError("the MCS search reported cancellation with no limit set")
    if not result.smartsString or result.numAtoms < 1:
        return None
    return Chem.MolFromSmarts(result.smartsString)


def charge_atom_ids(mol: Chem.Mol) -> Set[int]:
    return {int(a.GetIdx()) for a in mol.GetAtoms()
            if int(a.GetFormalCharge()) != 0}


def outside_components(mol: Chem.Mol, core: Set[int]) -> List[Set[int]]:
    """The connected pieces of what lies outside a candidate core."""
    seen: Set[int] = set()
    out: List[Set[int]] = []
    for start in (i for i in range(mol.GetNumAtoms()) if i not in core):
        if start in seen:
            continue
        component = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbour in mol.GetAtomWithIdx(current).GetNeighbors():
                index = int(neighbour.GetIdx())
                if index in core or index in component:
                    continue
                component.add(index)
                stack.append(index)
        seen |= component
        out.append(component)
    return out


def match_is_usable(mol: Chem.Mol, match: Sequence[int],
                    expected_charge: int) -> bool:
    """The three acceptance conditions, applied to one match.

    The charge test outside the core is taken per component, not per atom: a
    neutral component carrying a positive and a negative atom is still neutral,
    and excluding it would be stricter than the method claims to be.
    """
    core = set(map(int, match))
    charged = charge_atom_ids(mol)
    if not charged or not charged.issubset(core):
        return False
    if sum(int(mol.GetAtomWithIdx(i).GetFormalCharge())
           for i in core) != int(expected_charge):
        return False
    for component in outside_components(mol, core):
        if sum(int(mol.GetAtomWithIdx(i).GetFormalCharge())
               for i in component) != 0:
            return False
        if any(int(mol.GetAtomWithIdx(i).GetNumRadicalElectrons()) != 0
               for i in component):
            return False
    return True


def rank_matches(mol: Chem.Mol, query: Chem.Mol, expected_charge: int):
    """Decision 1 — which atoms of a member are its core.

    Fewest bonds leaving the match, then most heteroatoms, then most aromatic
    atoms, then the lowest atom indices. The last key exists only to make the
    choice unique and carries no chemical meaning.

    Symmetry-equivalent matches are collapsed before ranking, so a molecule
    with four equivalent substituents offers one match and not twenty-four.
    """
    matches = mol.GetSubstructMatches(query, uniquify=True, maxMatches=10000)
    usable = [tuple(map(int, m)) for m in matches
              if match_is_usable(mol, m, expected_charge)]
    if not usable:
        return None, None

    def score(match: Tuple[int, ...]):
        core = set(match)
        boundary = hetero = aromatic = 0
        for i in match:
            atom = mol.GetAtomWithIdx(i)
            hetero += int(atom.GetAtomicNum() not in (1, 6))
            aromatic += int(atom.GetIsAromatic())
            boundary += sum(int(nb.GetIdx()) not in core
                            for nb in atom.GetNeighbors())
        return (-boundary, hetero, aromatic, tuple(-x for x in match))

    ordered = sorted(usable, key=score, reverse=True)
    best = ordered[0]
    trace = _trace(ordered, score,
                   ("bonds_leaving_the_match", "heteroatoms", "aromatic_atoms",
                    "atom_index_order"))
    return best, trace


def _trace(ordered: Sequence[Any], keyfn,
           key_names: Sequence[str]) -> Dict[str, Any]:
    """Which key in a ranking discriminated, and how many were still tied."""
    if len(ordered) == 1:
        return {"key_that_discriminated": "only_one_option",
                "options_tied_at_that_key": 1, "options_considered": 1,
                "runner_up": ""}
    first, second = keyfn(ordered[0]), keyfn(ordered[1])
    which, tied = "all_keys_equal", len(ordered)
    for position, name in enumerate(key_names):
        if first[position] != second[position]:
            which = name
            tied = sum(1 for o in ordered
                       if keyfn(o)[:position] == first[:position])
            break
    return {"key_that_discriminated": which,
            "options_tied_at_that_key": tied,
            "options_considered": len(ordered),
            "runner_up": str(ordered[1])}


def build_exact_core_query(source: Chem.Mol,
                           match: Sequence[int]) -> Chem.Mol:
    """One query form: the member copied atom for atom.

    Built as a molecule and never round-tripped through SMARTS text, so the
    element, charge, aromaticity and bond type of every atom and bond survive
    exactly as the member carries them.
    """
    ids = tuple(map(int, match))
    inside = set(ids)
    builder = Chem.RWMol()
    position: Dict[int, int] = {}
    for old in ids:
        src = source.GetAtomWithIdx(old)
        atom = Chem.Atom(int(src.GetAtomicNum()))
        atom.SetFormalCharge(int(src.GetFormalCharge()))
        atom.SetIsAromatic(bool(src.GetIsAromatic()))
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        position[old] = int(builder.AddAtom(atom))
    for bond in source.GetBonds():
        i, j = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        if i in inside and j in inside:
            builder.AddBond(position[i], position[j], bond.GetBondType())
            added = builder.GetBondBetweenAtoms(position[i], position[j])
            if added is not None:
                added.SetIsAromatic(bool(bond.GetIsAromatic()))
                added.SetIsConjugated(bool(bond.GetIsConjugated()))
    out = builder.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def core_query_key(query: Chem.Mol) -> str:
    return Chem.MolToSmarts(query)


def big_enough(query: Chem.Mol) -> bool:
    return (query.GetNumAtoms() >= MIN_CORE_ATOMS
            and query.GetNumBonds() >= MIN_CORE_BONDS)


def has_ring(query: Chem.Mol) -> bool:
    probe = Chem.MolFromSmarts(Chem.MolToSmarts(query))
    if probe is None:
        return False
    Chem.GetSSSR(probe)
    return probe.GetRingInfo().NumRings() > 0


def captures_the_family_ring(query: Chem.Mol,
                             supported: Sequence[Chem.Mol]) -> bool:
    """A core must be cyclic if every ion it would hold is cyclic.

    The size floor stops a bare heteroatom winning on support, but it does not
    stop a three-atom acyclic fragment doing the same to a family whose every
    member carries a ring: the bicyclic cations, the cyclic sulfoniums and the
    cyclic phosphoniums each had their defining skeleton replaced by a short
    chain. A core that captures no ring where all its members have one is not
    the substructure those members share in any sense a chemist would accept.
    """
    if not supported:
        return True
    if any(m.GetRingInfo().NumRings() == 0 for m in supported):
        return True
    return has_ring(query)


@dataclass
class Evaluated:
    query: Chem.Mol = field(repr=False)
    support: Dict[int, Tuple[int, ...]]
    traces: Dict[int, Dict[str, Any]]
    form: str
    copied_from: int
    form_trace: Dict[str, Any]
    source_trace: Dict[str, Any]


def evaluate_candidate(raw: Chem.Mol,
                       records: Sequence[Dict[str, Any]],
                       mols: Dict[int, Chem.Mol]
                       ) -> Tuple[Optional[Evaluated], str]:
    """Decisions 2 and 3, and the support the winning form carries.

    Returns the winning form, or nothing and the name of the condition that
    declined every form. A single catch-all reason tells a reader that a
    candidate lost without telling them to what, and three different conditions
    can decline one here.

    Both query forms are kept and the one preserving the most support wins. The
    copied form can be strictly more specific than the shared motif, so a member
    differing only in a bond order or in aromaticity is silently dropped by it
    and a genuine common core is discarded as if it did not exist. Formal charge
    and charge localisation are enforced by the match test either way, so the
    flexible form never admits a neutral core.
    """
    provisional: Dict[int, Tuple[int, ...]] = {}
    provisional_traces: Dict[int, Dict[str, Any]] = {}
    for record in records:
        mol = mols[record["inventory_id"]]
        match, trace = rank_matches(mol, raw, record["required_charge"])
        if match is not None:
            provisional[record["inventory_id"]] = match
            provisional_traces[record["inventory_id"]] = trace
    if len(provisional) < MIN_SUPPORT:
        return None, ("the search pattern matches %d member(s) and a core "
                      "needs %d" % (len(provisional), MIN_SUPPORT))

    reached = [r for r in records if r["inventory_id"] in provisional]
    ordered = sorted(reached, key=lambda r: (r["heavy_atoms"],
                                             r["canonical_smiles"]))
    source = ordered[0]
    source_trace = _trace(ordered,
                          lambda r: (r["heavy_atoms"], r["canonical_smiles"]),
                          ("heavy_atoms", "canonical_smiles"))
    exact = build_exact_core_query(mols[source["inventory_id"]],
                                   provisional[source["inventory_id"]])

    def support_of(query: Chem.Mol):
        support: Dict[int, Tuple[int, ...]] = {}
        traces: Dict[int, Dict[str, Any]] = {}
        for record in records:
            mol = mols[record["inventory_id"]]
            match, trace = rank_matches(mol, query, record["required_charge"])
            if match is not None:
                support[record["inventory_id"]] = match
                traces[record["inventory_id"]] = trace
        return support, traces

    options: List[Tuple[Chem.Mol, Dict[int, Tuple[int, ...]], Dict[int, Dict[str, Any]], str]] = []
    built = 0
    if big_enough(exact):
        support, traces = support_of(exact)
        options.append((exact, support, traces, "member_copied_atom_for_atom"))
        built += 1
    if big_enough(raw):
        options.append((raw, dict(provisional), dict(provisional_traces),
                        "search_pattern"))
        built += 1
    if not built:
        return None, ("neither form reaches the floor of %d atoms and %d bonds"
                      % (MIN_CORE_ATOMS, MIN_CORE_BONDS))

    held = [o for o in options if len(o[1]) >= MIN_SUPPORT]
    if not held:
        return None, ("both forms clear the size floor and neither is shared "
                      "by %d members" % MIN_SUPPORT)
    cyclic = [o for o in held
              if captures_the_family_ring(o[0], [mols[i] for i in o[1]])]
    if not cyclic:
        return None, ("every member it holds carries a ring and no form of it "
                      "does, so it is not the substructure they share")
    options = cyclic

    def form_key(option):
        return (len(option[1]), option[0].GetNumAtoms(),
                option[0].GetNumBonds())

    options.sort(key=form_key, reverse=True)
    form_trace = _trace(options, form_key,
                        ("members_matched", "atoms", "bonds"))
    query, support, traces, form = options[0]
    return Evaluated(query=query, support=support, traces=traces, form=form,
                     copied_from=source["inventory_id"],
                     form_trace=form_trace, source_trace=source_trace), ""


def charge_environment(mol: Chem.Mol) -> str:
    """The local surroundings of every charged atom, as one grouping key."""
    parts = []
    for atom in mol.GetAtoms():
        if int(atom.GetFormalCharge()) == 0:
            continue
        neighbours = sorted((int(n.GetAtomicNum()), int(n.GetIsAromatic()),
                             int(n.GetFormalCharge()))
                            for n in atom.GetNeighbors())
        parts.append((int(atom.GetAtomicNum()), int(atom.GetFormalCharge()),
                      int(atom.GetIsAromatic()), int(atom.GetDegree()),
                      tuple(neighbours)))
    return repr(sorted(parts))
