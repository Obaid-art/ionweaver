
# -*- coding: utf-8 -*-
"""Stage 06 — the string controls, matched attempt for attempt.

For every graph candidate examined under the four atom-level
operators, one elementary edit on the same ion in each of SMILES and SELFIES.
Insert maps to inserting one token, Delete to deleting one, and Substitute and
Rebond each to replacing one. The attempt counts are therefore identical on all
three sides, which is the only thing that makes the comparison mean anything.

Both representations draw from the vocabulary of tokens the ion's own family
displays in that representation, so neither control is handed chemistry its
family never shows.

Silicon is declared before anything is scored. The default SELFIES
constraint table has no entry for it, so silicon falls through to the table's
catch-all of eight, and a four-valent element is given eight bonds' worth of
room. Ten of the inventory's elements are placeable and silicon is one of them,
so the control would otherwise be scored against a valence no chemistry has.

The parent-exclusion rule is the same one IonWeaver answers to. The previous
code counted a product identical to the parent as a success for the controls and
excluded it for IonWeaver, which flattered the controls on exactly the edits
most likely to return their input.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import selfies

from .common import STUDY_SCORING, RunError, Scoring, judge
from .operators import ATOM_LEVEL

#: Silicon's valence, declared because the default table does not carry it.
SILICON_VALENCE = 4

#: One token of a SMILES string. Two-letter elements come before the one-letter
#: organic subset so Cl is never read as C followed by l.
SMILES_TOKEN = re.compile(
    r"(\[[^\]]+\]|Br|Cl|Si|Se|se|@@|@|=|#|\\|/|%[0-9]{2}|[0-9]"
    r"|\(|\)|\.|-|\+|B|C|N|O|P|S|F|I|b|c|n|o|p|s|\*)")

#: Which elementary string edit stands opposite which operator. Substitute and
#: Rebond both map to a replacement: the string side offers three elementary
#: edits and the graph side has four atom-level operators, and the attempt
#: counts still match one to one because the mapping is per candidate.
EDIT_CLASS: Dict[str, str] = {
    "Insert": "insert",
    "Delete": "delete",
    "Substitute": "replace",
    "Rebond": "replace",
}

#: The mapping must name the atom-level operators and no others. Listing them
#: here and in `operators.ATOM_LEVEL` is two places to keep in step, so the
#: second one is checked against the first at import time instead.
if set(EDIT_CLASS) != set(ATOM_LEVEL):
    raise RunError(
        "the string controls map %s and the atom-level operators are %s. The "
        "matched comparison runs on the atom-level four, so every one of them "
        "needs an elementary string edit opposite it and nothing else may have "
        "one." % (sorted(EDIT_CLASS), sorted(ATOM_LEVEL)))

CONTROL_COLUMNS = (
    "attempt_id", "method", "inventory_id", "family_id",
    "matched_examination_id", "matched_operator", "edit_class", "position",
    "token_before", "token_after", "parent_string", "edited_string",
    "product_smiles", "sanitises", "connected", "closed_shell",
    "required_charge", "core_retained", "charge_on_core",
    "identical_to_parent", "accepted", "reason")

_DECLARED = False


def declare_silicon() -> Dict[str, int]:
    """Put silicon in the constraint table, once, before anything is scored."""
    global _DECLARED
    table = dict(selfies.get_semantic_constraints())
    if table.get("Si") != SILICON_VALENCE:
        table["Si"] = SILICON_VALENCE
        selfies.set_semantic_constraints(table)
    _DECLARED = True
    return dict(selfies.get_semantic_constraints())


def smiles_tokens(text: str) -> List[str]:
    """Split a SMILES string into tokens, and prove the split is lossless."""
    found = SMILES_TOKEN.findall(text)
    if "".join(found) != text:
        raise RunError(
            "the SMILES tokeniser did not reproduce %r from its tokens; it "
            "produced %r. A control that edits a token it cannot write back is "
            "not editing the molecule." % (text, "".join(found)))
    return found


def selfies_tokens(text: str) -> List[str]:
    return list(selfies.split_selfies(text))


def vocabularies(ions: Sequence[Dict[str, Any]]
                 ) -> Dict[str, Dict[str, Tuple[str, ...]]]:
    """Per family, the token vocabulary each representation displays."""
    if not _DECLARED:
        declare_silicon()
    smiles_vocab: Dict[str, Set[str]] = {}
    selfies_vocab: Dict[str, Set[str]] = {}
    for record in ions:
        family = record["family_id"]
        text = record["canonical_smiles"]
        smiles_vocab.setdefault(family, set()).update(smiles_tokens(text))
        encoded = record.get("selfies") or ""
        if not encoded:
            try:
                encoded = selfies.encoder(text)
            except Exception:                              # noqa: BLE001
                encoded = ""
        if encoded:
            selfies_vocab.setdefault(family, set()).update(
                selfies_tokens(encoded))
    return {"smiles": {k: tuple(sorted(v)) for k, v in smiles_vocab.items()},
            "selfies": {k: tuple(sorted(v))
                        for k, v in selfies_vocab.items()}}


@dataclass(frozen=True)
class StringEdit:
    edit_class: str
    position: int
    token_before: str
    token_after: str

    def apply(self, tokens: Sequence[str]) -> List[str]:
        out = list(tokens)
        if self.edit_class == "insert":
            out.insert(self.position, self.token_after)
        elif self.edit_class == "delete":
            del out[self.position]
        elif self.edit_class == "replace":
            out[self.position] = self.token_after
        else:
            raise RunError("unknown elementary string edit %r"
                           % self.edit_class)
        return out


def options(tokens: Sequence[str], vocabulary: Sequence[str],
            edit_class: str) -> List[StringEdit]:
    """Every elementary edit of this class available on this string."""
    out: List[StringEdit] = []
    if edit_class == "insert":
        for position in range(len(tokens) + 1):
            for token in vocabulary:
                out.append(StringEdit("insert", position, "", token))
    elif edit_class == "delete":
        for position in range(len(tokens)):
            out.append(StringEdit("delete", position, tokens[position], ""))
    elif edit_class == "replace":
        for position in range(len(tokens)):
            for token in vocabulary:
                if token == tokens[position]:
                    continue
                out.append(StringEdit("replace", position, tokens[position],
                                      token))
    else:
        raise RunError("unknown elementary string edit %r" % edit_class)
    return out


def shuffle_key(seed: int, inventory_id: int, method: str, edit_class: str,
                edit: StringEdit) -> str:
    """The same keyed shuffle the graph side uses, so neither is ordered by an
    arithmetic progression."""
    return hashlib.blake2b(
        ("%d|%d|%s|%s|%d|%s|%s"
         % (seed, inventory_id, method, edit_class, edit.position,
            edit.token_before, edit.token_after)).encode("utf-8"),
        digest_size=16).hexdigest()


def attempts_for(inventory_id: int, family_id: str, method: str,
                 parent_string: str, tokens: Sequence[str],
                 vocabulary: Sequence[str], wanted: Sequence[Tuple[int, str]],
                 seed: int) -> List[Tuple[int, str, StringEdit]]:
    """One string edit per graph candidate, of the matching class.

    `wanted` is the ordered list of (examination id, operator) the graph side
    examined under the atom-level four. Each gets one string edit of the class
    that operator maps to, drawn without replacement from the options of that
    class in keyed order.
    """
    ordered: Dict[str, List[StringEdit]] = {}
    used: Dict[str, int] = {}
    out: List[Tuple[int, str, StringEdit]] = []
    for examination_id, operator in wanted:
        edit_class = EDIT_CLASS[operator]
        if edit_class not in ordered:
            available = options(tokens, vocabulary, edit_class)
            available.sort(key=lambda e: shuffle_key(seed, inventory_id,
                                                     method, edit_class, e))
            ordered[edit_class] = available
            used[edit_class] = 0
        index = used[edit_class]
        if index >= len(ordered[edit_class]):
            continue           # this string offers fewer; the count is recorded
        used[edit_class] = index + 1
        out.append((examination_id, operator, ordered[edit_class][index]))
    del parent_string, family_id
    return out


def score(method: str, inventory_id: int, family_id: str,
          parent_smiles: str, parent_string: str, tokens: Sequence[str],
          vocabulary: Sequence[str], wanted: Sequence[Tuple[int, str]],
          core_pattern, required_charge: int, seed: int,
          scoring: Scoring = STUDY_SCORING) -> List[Dict[str, Any]]:
    """Run and judge one ion's control attempts in one representation."""
    if not _DECLARED:
        declare_silicon()
    rows: List[Dict[str, Any]] = []
    for examination_id, operator, edit in attempts_for(
            inventory_id, family_id, method, parent_string, tokens,
            vocabulary, wanted, seed):
        edited = "".join(edit.apply(tokens))
        if method == "selfies":
            try:
                product = selfies.decoder(edited)
            except Exception:                              # noqa: BLE001
                product = ""
        else:
            product = edited
        verdict = judge(product, parent_smiles, core_pattern, required_charge,
                        scoring) if product else None
        row = {
            "method": method, "inventory_id": inventory_id,
            "family_id": family_id,
            "matched_examination_id": examination_id,
            "matched_operator": operator, "edit_class": edit.edit_class,
            "position": edit.position, "token_before": edit.token_before,
            "token_after": edit.token_after,
            "parent_string": parent_string, "edited_string": edited,
            "product_smiles": product}
        if verdict is None:
            row.update({"sanitises": 0, "connected": 0, "closed_shell": 0,
                        "required_charge": 0, "core_retained": 0,
                        "charge_on_core": 0, "identical_to_parent": 0,
                        "accepted": 0,
                        "reason": "the edited string does not decode"})
        else:
            row.update(verdict.as_row())
        rows.append(row)
    return rows
