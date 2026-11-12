"""The published state, rebuilt from the archive tables at start-up.

The submitted run wrote every stage to a table.  Three of those stages are
cheap to recompute and one is not, so the server:

* re-curates the 2,431 source entries with the engine's own ``curation.curate``
  (a few seconds) and checks that it lands on the same 1,241 ions with the same
  identifiers the archive carries;
* re-learns the atom, bond and ring environment registries from those ions with
  ``registry.build`` (under a second);
* reads the 44 cores and 1,222 core assignments from the archive instead of
  re-running the maximum-common-substructure search (35 s in the run), then
  rebuilds the editable set and the fragment library with ``fragments``;
* checks the rebuilt editable set, the library and the registries against the
  archive's own tables, so a server that does not reproduce the run refuses to
  start rather than serve something else.
"""
from __future__ import annotations

import collections
import csv
import functools
import json
import time
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional

from rdkit import Chem

from . import DATA_DIR
from iw import curation, fragments, registry

__all__ = ["Published", "load_published", "CoreModel", "SEED"]

#: The run seed of the submitted study, used so a reference ion's candidates
#: come out in the published order.
SEED = 20261005


def _read(name: str) -> List[Dict[str, str]]:
    with open(DATA_DIR / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _json(name: str) -> Any:
    return json.loads((DATA_DIR / name).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class CoreModel:
    family_id: str
    family_label: str
    ion_type: str
    grammar_id: str
    pattern: str
    atoms: int
    bonds: int
    charged_atoms: int
    support: int
    support_fraction: float
    copied_from: int
    query: Any = field(compare=False, hash=False, repr=False)

    def as_dict(self) -> Dict[str, Any]:
        return {"family_id": self.family_id, "family_label": self.family_label,
                "ion_type": self.ion_type, "grammar_id": self.grammar_id,
                "pattern": self.pattern, "atoms": self.atoms, "bonds": self.bonds,
                "charged_atoms": self.charged_atoms, "support": self.support,
                "support_fraction": self.support_fraction, "copied_from": self.copied_from}


class Published:
    """Everything the engine needs, plus the archive's own numbers."""

    def __init__(self) -> None:
        t0 = time.perf_counter()
        self.checks: Dict[str, Any] = {}

        # -- 02 curation, recomputed from the raw columns of the audit ---------
        audit_rows = _read("02_curation__curation_audit.csv")
        entries = [{
            "raw_index": int(r["raw_index"]), "source_workbook": r["source_workbook"],
            "source_sheet": r["source_sheet"], "source_row": int(r["source_row"]),
            "ion_type_raw": r["ion_type_raw"], "abbreviation_raw": r["abbreviation_raw"],
            "chemical_name_raw": r["chemical_name_raw"], "smiles_raw": r["smiles_raw"],
            "declared_charge_raw": r["declared_charge_raw"], "family_raw": r["family_raw"],
        } for r in audit_rows]
        self.audit, self.inventory = curation.curate(entries)
        self.curation_summary = curation.summary(self.audit, self.inventory)
        published_retained = [(int(r["raw_index"]), r["canonical_smiles"]) for r in audit_rows
                              if r["status"] == "retained"]
        rebuilt = [(r["raw_index"], r["canonical_smiles"]) for r in self.inventory]
        self.checks["curation_reproduced"] = published_retained == rebuilt
        self.by_inventory_id: Dict[int, Dict[str, Any]] = {r["inventory_id"]: r for r in self.inventory}

        # -- 01 registries ---------------------------------------------------------
        reg = registry.build(self.inventory)
        self.observed: FrozenSet = frozenset(reg["observed"])
        self.observed_bonds: FrozenSet = frozenset(reg["observed_bonds"])
        self.creatable: FrozenSet = frozenset(reg["creatable"])
        self.elements: List[int] = sorted({e[0] for e in reg["creatable"]})
        self.registry_rows = {"atoms": reg["atom_rows"], "bonds": reg["bond_rows"], "rings": reg["ring_rows"]}
        self.checks["atom_environments_reproduced"] = (
            len(reg["atom_rows"]) == len(_read("01_inputs__observed_atom_environments.csv")))
        self.checks["bond_environments_reproduced"] = (
            len(reg["bond_rows"]) == len(_read("01_inputs__observed_bond_environments.csv")))

        # -- 03 cores, read from the archive ----------------------------------------
        labels = {r["family_id"]: r["family_label"] for r in self.inventory}
        models = []
        for r in _read("03_cores__family_core_models.csv"):
            models.append({
                "family_id": r["family_id"], "grammar_id": r["grammar_id"], "layer": int(r["layer"]),
                "pattern": r["pattern"], "atoms": int(r["atoms"]), "bonds": int(r["bonds"]),
                "charged_atoms": int(r["charged_atoms"]), "support": int(r["support"]),
                "support_fraction_of_family": float(r["support_fraction_of_family"]),
                "copied_from_inventory_id": int(r["copied_from_inventory_id"]),
                "query_form": r["query_form"]})
        assignments = [{
            "inventory_id": int(r["inventory_id"]), "family_id": r["family_id"],
            "grammar_id": r["grammar_id"], "pattern": r["pattern"],
            "matched_atoms": r["matched_atoms"],
            "usable_matches_available": int(r["usable_matches_available"]),
            "match_chosen_by_key": r["match_chosen_by_key"],
        } for r in _read("03_cores__source_core_assignments.csv")]
        self.cores = {"models": models, "assignments": assignments}
        self.core_models: List[CoreModel] = [CoreModel(
            family_id=m["family_id"], family_label=labels.get(m["family_id"], m["family_id"]),
            ion_type="cation" if m["family_id"].startswith("cation") else "anion",
            grammar_id=m["grammar_id"], pattern=m["pattern"], atoms=m["atoms"], bonds=m["bonds"],
            charged_atoms=m["charged_atoms"], support=m["support"],
            support_fraction=m["support_fraction_of_family"], copied_from=m["copied_from_inventory_id"],
            query=Chem.MolFromSmarts(m["pattern"])) for m in models]
        self.uncovered = _read("03_cores__ions_without_a_core.csv")
        self.coverage = _read("03_cores__family_core_coverage.csv")

        # -- 04 the editable set and the fragment library ---------------------------
        self.starting: List[Dict[str, Any]] = fragments.starting_ions(self.cores, self.inventory)
        self.by_id: Dict[int, Dict[str, Any]] = {r["inventory_id"]: r for r in self.starting}
        archived = _read("03_cores__starting_ions.csv")
        self.checks["starting_ions_reproduced"] = (
            [(int(a["inventory_id"]), a["core_atoms"], a["core_pattern"], int(a["body_heavy_atoms"]))
             for a in archived]
            == [(r["inventory_id"], r["core_atoms"], r["core_pattern"], r["body_heavy_atoms"])
                for r in self.starting])
        built = fragments.build(self.starting)
        self.library: List[Dict[str, Any]] = built["library"]
        self.fragments: List[Dict[str, Any]] = built["fragments"]
        self.instances: List[Dict[str, Any]] = built["instances"]
        self.library_summary = fragments.summary(built)
        self.checks["library_reproduced"] = (
            [r["component_key"] for r in _read("04_fragments__component_library.csv")]
            == [r["component_key"] for r in self.library])
        self.checks["fragments_reproduced"] = (
            [r["fragment_key"] for r in _read("04_fragments__terminal_fragments.csv")]
            == [r["fragment_key"] for r in self.fragments])
        self.silent = _read("03_cores__ions_without_an_editable_site.csv")

        # -- the archive's own numbers, for the page ----------------------------------
        self.headline = _json("10_aggregates__headline_results.json")
        self.sensitivity = _read("05_single_edits__budget_sensitivity.csv")
        self.string_controls = _read("06_string_controls__string_control_summary.csv")
        self.per_ion_operator = _read("05_single_edits__per_ion_operator_summary.csv")
        self.smiles_set = {r["canonical_smiles"] for r in self.inventory}
        self.load_seconds = round(time.perf_counter() - t0, 2)
        self.checks["all"] = all(v for k, v in self.checks.items() if k != "all")

    # -- lookups -------------------------------------------------------------------

    def record(self, inventory_id: int) -> Dict[str, Any]:
        try:
            return self.by_id[int(inventory_id)]
        except KeyError as exc:
            raise KeyError(f"no editable reference ion with id {inventory_id}") from exc

    def record_by_smiles(self, canonical: str) -> Optional[Dict[str, Any]]:
        return next((r for r in self.starting if r["canonical_smiles"] == canonical), None)

    @functools.cached_property
    def families(self) -> List[Dict[str, Any]]:
        """One row per family: label, type, members, editable members, cores."""
        rows: Dict[str, Dict[str, Any]] = {}
        members = collections.Counter(r["family_id"] for r in self.inventory)
        for r in self.inventory:
            rows.setdefault(r["family_id"], {
                "family_id": r["family_id"], "family_label": r["family_label"],
                "ion_type": r["ion_type"], "members": members[r["family_id"]],
                "editable": 0, "cores": [], "representative_id": None, "_rep_body": 10 ** 9})
        for s in self.starting:
            row = rows[s["family_id"]]
            row["editable"] += 1
            if s["grammar_id"] == "G0" and 0 < s["body_heavy_atoms"] < row["_rep_body"]:
                row["_rep_body"] = s["body_heavy_atoms"]
                row["representative_id"] = s["inventory_id"]
        for m in self.core_models:
            rows[m.family_id]["cores"].append(m.as_dict())
        for row in rows.values():
            row.pop("_rep_body")
            if row["representative_id"] is None:
                first = next((s for s in self.starting if s["family_id"] == row["family_id"]), None)
                row["representative_id"] = first["inventory_id"] if first else None
            row["cores"].sort(key=lambda c: c["grammar_id"])
        return sorted(rows.values(), key=lambda r: (r["ion_type"], r["family_label"]))

    def family(self, family_id: str) -> Optional[Dict[str, Any]]:
        return next((f for f in self.families if f["family_id"] == family_id), None)

    def ions_of(self, family_id: str) -> List[Dict[str, Any]]:
        return sorted((s for s in self.starting if s["family_id"] == family_id),
                      key=lambda s: (s["grammar_id"], s["heavy_atoms"], s["inventory_id"]))

    def counts(self) -> Dict[str, int]:
        return {"source_entries": len(self.audit), "ions": len(self.inventory),
                "families": len(self.families), "cores": len(self.core_models),
                "editable_ions": len(self.starting), "library_components": len(self.library),
                "terminal_fragments": len(self.fragments),
                "atom_environments": len(self.registry_rows["atoms"]),
                "bond_environments": len(self.registry_rows["bonds"])}


@functools.lru_cache(maxsize=1)
def load_published() -> Published:
    return Published()
