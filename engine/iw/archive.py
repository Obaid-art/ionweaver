
# -*- coding: utf-8 -*-
"""The archive. Every stage's records on disk, in the layout spec E2 names.

Nothing may exist only in a summary, so each stage writes its rows before any
count is computed from them, and every count a later stage reports is read back
from these files. An empty table is written with its header, because an empty
file is a result.
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Sequence

from .common import write_csv, write_json

#: the stage directories, in the order the run fills them
STAGES = ("01_inputs", "02_curation", "03_cores", "04_fragments",
          "05_single_edits", "06_string_controls", "07_sequential",
          "08_comparison", "09_validation", "10_aggregates")

CURATION_COLUMNS = (
    "raw_index", "source_workbook", "source_sheet", "source_row",
    "ion_type_raw", "abbreviation_raw", "chemical_name_raw", "smiles_raw",
    "declared_charge_raw", "family_raw", "family_id", "required_charge",
    "in_scope", "scope_reason", "structure_check", "status",
    "canonical_smiles", "heavy_atoms")

#: Imported, never copied. A second copy of a column list is how the hydrogen
#: count came out as an empty column the first time this stage ran.
from .registry import (ATOM_ENV_COLUMNS, BOND_ENV_COLUMNS,  # noqa: E402
                       RING_ENV_COLUMNS)

CORE_MODEL_COLUMNS = (
    "family_id", "grammar_id", "layer", "pattern", "atoms", "bonds",
    "charged_atoms", "support", "support_fraction_of_family",
    "copied_from_inventory_id", "query_form")

CORE_DIAGNOSTIC_COLUMNS = (
    "family_id", "layer", "pass", "support_threshold", "search_pattern",
    "pattern", "atoms", "bonds", "charged_atoms", "support", "form",
    "copied_from", "status", "reason")

CORE_DECISION_COLUMNS = (
    "family_id", "layer", "grammar_id", "decision_point", "rule", "chosen",
    "runner_up", "key_that_discriminated", "options_tied_at_that_key",
    "options_considered")

ASSIGNMENT_COLUMNS = (
    "inventory_id", "family_id", "grammar_id", "pattern", "matched_atoms",
    "usable_matches_available", "match_chosen_by_key")

COVERAGE_COLUMNS = ("family_id", "members", "members_assigned", "cores",
                    "coverage")

STARTING_COLUMNS = (
    "inventory_id", "family_id", "family_label", "canonical_smiles",
    "grammar_id", "core_pattern", "core_atoms", "core_atom_count",
    "ring_shared_atoms", "ring_shared_count", "heavy_atoms",
    "body_heavy_atoms", "body_components")

UNCOVERED_COLUMNS = ("inventory_id", "family_id", "canonical_smiles", "reason")

RING_SHARED_COLUMNS = (
    "inventory_id", "family_id", "canonical_smiles", "core_pattern",
    "core_atoms", "ring_shared_atoms", "ring_shared_count", "heavy_atoms",
    "reason")


def prepare(root: str) -> None:
    """Make the stage directories. Writing is per stage, not per run."""
    for stage in STAGES:
        os.makedirs(os.path.join(root, stage), exist_ok=True)


def engine_sources() -> List[Dict[str, Any]]:
    """The engine's own source, hashed as it ran.

    Taken from each loaded module's loader rather than from a directory on
    disk, because the engine ships two ways: as files, and carried inside a
    single Colab script. Hashing the loader's source covers both, and it
    hashes what actually ran instead of what happens to be lying beside the
    archive. A module whose source cannot be recovered is listed as such; it is
    not silently left out of the manifest.
    """
    import hashlib
    import sys

    out: List[Dict[str, Any]] = []
    for name in sorted(sys.modules):
        # the package itself is an empty namespace with no code of its own
        if not name.startswith("iw."):
            continue
        module = sys.modules[name]
        loader = getattr(module, "__loader__", None)
        source = None
        if loader is not None and hasattr(loader, "get_source"):
            try:
                source = loader.get_source(name)
            except Exception:                              # noqa: BLE001
                source = None
        if source is None:
            path = getattr(module, "__file__", None)
            if path and os.path.exists(path):
                source = open(path, encoding="utf-8").read()
        if source is None:
            out.append({"file": "engine/%s" % name, "bytes": "",
                        "sha256": "source not recoverable from this module"})
            continue
        raw = source.encode("utf-8")
        out.append({"file": "engine/%s" % name, "bytes": len(raw),
                    "sha256": hashlib.sha256(raw).hexdigest()})
    return out


def at(root: str, stage: str, name: str) -> str:
    return os.path.join(root, stage, name)


def write_inputs(root: str, provenance: Dict[str, Any],
                 registry: Dict[str, Any]) -> Dict[str, int]:
    """01_inputs."""
    write_json(at(root, "01_inputs", "input_provenance.json"), provenance)
    return {
        "observed_atom_environments.csv": write_csv(
            at(root, "01_inputs", "observed_atom_environments.csv"),
            registry["atom_rows"], ATOM_ENV_COLUMNS),
        "observed_bond_environments.csv": write_csv(
            at(root, "01_inputs", "observed_bond_environments.csv"),
            registry["bond_rows"], BOND_ENV_COLUMNS),
        "observed_ring_environments.csv": write_csv(
            at(root, "01_inputs", "observed_ring_environments.csv"),
            registry["ring_rows"], RING_ENV_COLUMNS),
    }


def write_curation(root: str, audit: Sequence[Dict[str, Any]],
                   summary: Dict[str, Any]) -> Dict[str, int]:
    """02_curation."""
    write_json(at(root, "02_curation", "curation_summary.json"), summary)
    return {"curation_audit.csv": write_csv(
        at(root, "02_curation", "curation_audit.csv"), audit,
        CURATION_COLUMNS)}


def coverage_rows(inventory: Sequence[Dict[str, Any]],
                  models: Sequence[Dict[str, Any]],
                  assignments: Sequence[Dict[str, Any]]
                  ) -> List[Dict[str, Any]]:
    """Per family: members, how many a core reached, and how many cores."""
    import collections
    members = collections.Counter(r["family_id"] for r in inventory)
    assigned = collections.Counter(r["family_id"] for r in assignments)
    cores = collections.Counter(m["family_id"] for m in models)
    return [{"family_id": family,
             "members": members[family],
             "members_assigned": assigned[family],
             "cores": cores[family],
             "coverage": round(assigned[family] / members[family], 4)}
            for family in sorted(members)]


def write_cores(root: str, cores: Dict[str, Any],
                inventory: Sequence[Dict[str, Any]],
                starting: Sequence[Dict[str, Any]]) -> Dict[str, int]:
    """03_cores."""
    write_json(at(root, "03_cores", "core_settings.json"), cores["settings"])
    floor = support_floor(cores["models"], inventory,
                          cores["settings"]["min_support"])
    written = {
        "family_core_models.csv": write_csv(
            at(root, "03_cores", "family_core_models.csv"),
            cores["models"], CORE_MODEL_COLUMNS),
        "core_candidate_diagnostics.csv": write_csv(
            at(root, "03_cores", "core_candidate_diagnostics.csv"),
            cores["diagnostics"], CORE_DIAGNOSTIC_COLUMNS),
        "core_selection_decisions.csv": write_csv(
            at(root, "03_cores", "core_selection_decisions.csv"),
            cores["decisions"], CORE_DECISION_COLUMNS),
        "source_core_assignments.csv": write_csv(
            at(root, "03_cores", "source_core_assignments.csv"),
            cores["assignments"], ASSIGNMENT_COLUMNS),
        "family_core_coverage.csv": write_csv(
            at(root, "03_cores", "family_core_coverage.csv"),
            coverage_rows(inventory, cores["models"], cores["assignments"]),
            COVERAGE_COLUMNS),
        "ions_without_a_core.csv": write_csv(
            at(root, "03_cores", "ions_without_a_core.csv"),
            cores["uncovered"], UNCOVERED_COLUMNS),
        "starting_ions.csv": write_csv(
            at(root, "03_cores", "starting_ions.csv"), starting,
            STARTING_COLUMNS),
        "core_support_distribution.csv": write_csv(
            at(root, "03_cores", "core_support_distribution.csv"),
            floor["distribution"], SUPPORT_DISTRIBUTION_COLUMNS),
        "support_floor_cost.csv": write_csv(
            at(root, "03_cores", "support_floor_cost.csv"), floor["costs"],
            SUPPORT_FLOOR_COLUMNS),
    }
    write_json(at(root, "03_cores", "support_floor_finding.json"),
               floor["finding"])
    # `ions_without_an_editable_site.csv` is not written here. An ion with an
    # empty body still has editable positions, because a graft or an appended
    # insertion lands at a free valence on a core atom or on an atom held in a
    # ring the core covers in part. Only the enumeration knows which ions no
    # operator can touch, so only the stage that enumerates writes that file.
    return written


def write_fragments(root: str, result: Dict[str, Any],
                    summary: Dict[str, Any]) -> Dict[str, int]:
    """04_fragments."""
    from .fragments import (FRAGMENT_COLUMNS, INSTANCE_COLUMNS,
                            LIBRARY_COLUMNS)
    write_json(at(root, "04_fragments", "library_settings.json"),
               result["settings"])
    write_json(at(root, "04_fragments", "library_summary.json"), summary)
    return {
        "component_instances.csv": write_csv(
            at(root, "04_fragments", "component_instances.csv"),
            result["instances"], INSTANCE_COLUMNS),
        "component_library.csv": write_csv(
            at(root, "04_fragments", "component_library.csv"),
            result["library"], LIBRARY_COLUMNS),
        "terminal_fragments.csv": write_csv(
            at(root, "04_fragments", "terminal_fragments.csv"),
            result["fragments"], FRAGMENT_COLUMNS),
        "ring_shared_atoms.csv": write_csv(
            at(root, "04_fragments", "ring_shared_atoms.csv"),
            result["ring_shared"], RING_SHARED_COLUMNS),
    }


SUPPORT_DISTRIBUTION_COLUMNS = ("support", "cores", "cores_at_or_below",
                                "ions_held_at_this_support")

SUPPORT_FLOOR_COLUMNS = ("floor", "cores_dropped", "ions_losing_their_core",
                         "families_left_with_no_core_at_all",
                         "which_families", "their_member_counts")


def support_floor(models: Sequence[Dict[str, Any]],
                  inventory: Sequence[Dict[str, Any]],
                  floor_in_force: int) -> Dict[str, Any]:
    """How much support each core has, and what a higher floor would cost.

    A core supported by two members is a pairwise overlap, and a reader is
    entitled to ask why two is allowed. The answer has to come from the data,
    so the run reports the support distribution and what every candidate floor
    would remove, including the families it would leave with no core at all.
    """
    import collections

    members = collections.Counter(r["family_id"] for r in inventory)
    distribution = collections.Counter(m["support"] for m in models)

    rows: List[Dict[str, Any]] = []
    running = 0
    for support in sorted(distribution):
        running += distribution[support]
        rows.append({"support": support, "cores": distribution[support],
                     "cores_at_or_below": running,
                     "ions_held_at_this_support":
                         support * distribution[support]})

    costs: List[Dict[str, Any]] = []
    for floor in range(2, 6):
        dropped = [m for m in models if m["support"] < floor]
        kept = {m["family_id"] for m in models if m["support"] >= floor}
        emptied = sorted({m["family_id"] for m in dropped} - kept)
        costs.append({
            "floor": floor, "cores_dropped": len(dropped),
            "ions_losing_their_core": sum(m["support"] for m in dropped),
            "families_left_with_no_core_at_all": len(emptied),
            "which_families": " ".join(emptied),
            "their_member_counts": " ".join(str(members[f]) for f in emptied)})

    supports = sorted(distribution)
    gaps = sorted(((supports[i + 1] - supports[i], supports[i],
                    supports[i + 1]) for i in range(len(supports) - 1)),
                  reverse=True)
    lowest = [m for m in models if m["support"] == floor_in_force]
    return {
        "distribution": rows, "costs": costs,
        "finding": {
            "floor_in_force": floor_in_force,
            "supports_observed": supports,
            "largest_gaps": [{"gap": g, "between": a, "and": b}
                             for g, a, b in gaps[:5]],
            "cores_at_the_floor": [
                {"family_id": m["family_id"],
                 "grammar_id": m["grammar_id"], "pattern": m["pattern"],
                 "support": m["support"]} for m in lowest],
            "rule": "the floor is the fewest members a core may be supported "
                    "by. It is set from this distribution: a floor placed "
                    "where the distribution has no gap is an arbitrary cut, "
                    "and the cost table says which families each candidate "
                    "floor would leave with no core at all.",
        },
    }
