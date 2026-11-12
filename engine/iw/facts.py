
# -*- coding: utf-8 -*-
"""Stage 10 — the fact table, read back out of the archive.

Every quantity the paper may quote is computed here by reading the
written files, not by carrying a number forward from the stage that produced
it. That is the whole point: a number cannot reach the manuscript without
first having been written to a file and then read out of it, so a figure in the
paper and a figure in the deposit cannot drift apart.

Each fact carries the file it came from and the rule that produced it, in words,
so a reader can recompute it. A fact whose file is missing is reported as
missing and not as zero.
"""
from __future__ import annotations

import collections
import csv
import json
import os
import statistics
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

Rows = List[Dict[str, str]]


def read(root: str, relative: str) -> Optional[Rows]:
    path = os.path.join(root, relative)
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _int(row: Dict[str, str], key: str) -> int:
    value = row.get(key, "")
    return int(value) if value not in ("", None) else 0


def _float(row: Dict[str, str], key: str) -> Optional[float]:
    value = row.get(key, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


#: (name, file, rule, how). `how` is handed the file's rows and returns a
#: value. Nothing here reaches into a pickle or a stage's return value.
TABLE: Sequence[Tuple[str, str, str, Callable[[Rows], Any]]] = (
    ("source_entries_read",
     "02_curation/curation_audit.csv",
     "every row of the file, one per source entry",
     len),
    ("retained_ions",
     "02_curation/curation_audit.csv",
     "rows whose status is retained",
     lambda rows: sum(1 for r in rows if r["status"] == "retained")),
    ("families_retained",
     "02_curation/curation_audit.csv",
     "distinct family_id among the retained rows",
     lambda rows: len({r["family_id"] for r in rows
                       if r["status"] == "retained"})),
    ("entries_excluded_as_duplicates",
     "02_curation/curation_audit.csv",
     "rows whose status is excluded_duplicate_in_family",
     lambda rows: sum(1 for r in rows
                      if r["status"] == "excluded_duplicate_in_family")),
    ("entries_excluded_out_of_scope",
     "02_curation/curation_audit.csv",
     "rows whose status is excluded_out_of_scope",
     lambda rows: sum(1 for r in rows
                      if r["status"] == "excluded_out_of_scope")),
    ("entries_excluded_on_declared_charge",
     "02_curation/curation_audit.csv",
     "rows whose status is excluded_declared_charge_not_required",
     lambda rows: sum(
         1 for r in rows
         if r["status"] == "excluded_declared_charge_not_required")),

    ("atom_environments_observed",
     "01_inputs/observed_atom_environments.csv",
     "every row of the file, one per environment the inventory displays",
     len),
    ("bond_environments_observed",
     "01_inputs/observed_bond_environments.csv",
     "every row of the file, one per bonded pair the inventory displays",
     lambda rows: len(rows)),
    ("atom_environments_creatable",
     "01_inputs/observed_atom_environments.csv",
     "rows whose creatable flag is 1",
     lambda rows: sum(1 for r in rows if r["creatable"] == "1")),
    ("elements_an_edit_may_place",
     "01_inputs/observed_atom_environments.csv",
     "distinct element among the creatable rows",
     lambda rows: len({r["element"] for r in rows if r["creatable"] == "1"})),
    ("ring_environments_observed",
     "01_inputs/observed_ring_environments.csv",
     "every row of the file, one per (bond order, ring size) pair",
     len),

    ("cores_accepted",
     "03_cores/family_core_models.csv",
     "every row of the file, one per accepted core",
     len),
    ("families_with_more_than_one_core",
     "03_cores/family_core_coverage.csv",
     "rows whose core count exceeds one",
     lambda rows: sum(1 for r in rows if _int(r, "cores") > 1)),
    ("ions_a_core_reached",
     "03_cores/source_core_assignments.csv",
     "every row of the file, one per ion a core matched",
     len),
    ("ions_without_a_core",
     "03_cores/ions_without_a_core.csv",
     "every row of the file",
     len),
    ("ions_with_a_core_and_no_site_for_any_operator",
     "03_cores/ions_without_an_editable_site.csv",
     "every row of the file",
     len),
    ("core_selection_decisions_recorded",
     "03_cores/core_selection_decisions.csv",
     "every row of the file, one per decision point per accepted core",
     len),
    ("core_candidates_considered",
     "03_cores/core_candidate_diagnostics.csv",
     "every row of the file",
     len),
    ("cores_supported_by_the_fewest_members_allowed",
     "03_cores/core_support_distribution.csv",
     "the cores count on the row whose support is the floor in force",
     lambda rows: min([_int(r, "cores") for r in rows
                       if _int(r, "support") == min(_int(x, "support")
                                                    for x in rows)]
                      or [0])),
    ("families_a_floor_of_three_would_leave_without_any_core",
     "03_cores/support_floor_cost.csv",
     "the families_left_with_no_core_at_all on the floor-of-three row",
     lambda rows: next((_int(r, "families_left_with_no_core_at_all")
                        for r in rows if _int(r, "floor") == 3), None)),

    ("component_instances",
     "04_fragments/component_instances.csv",
     "every row of the file, one per body component of one ion",
     len),
    ("library_components",
     "04_fragments/component_library.csv",
     "every row of the file",
     len),
    ("library_components_of_one_heavy_atom",
     "04_fragments/component_library.csv",
     "rows whose heavy_atoms is one",
     lambda rows: sum(1 for r in rows if _int(r, "heavy_atoms") == 1)),
    ("library_components_attached_at_more_than_one_point",
     "04_fragments/component_library.csv",
     "rows whose attachment_count exceeds one",
     lambda rows: sum(1 for r in rows if _int(r, "attachment_count") > 1)),
    ("largest_library_component_heavy_atoms",
     "04_fragments/component_library.csv",
     "the largest heavy_atoms in the file",
     lambda rows: max([_int(r, "heavy_atoms") for r in rows] or [0])),
    ("terminal_fragments_at_two_or_more_heavy_atoms",
     "04_fragments/terminal_fragments.csv",
     "every row of the file",
     len),
    ("terminal_fragments_that_are_not_library_components",
     "04_fragments/terminal_fragments.csv",
     "rows whose is_itself_a_library_component flag is 0",
     lambda rows: sum(1 for r in rows
                      if r["is_itself_a_library_component"] == "0")),
    ("ions_holding_a_ring_the_core_covers_in_part",
     "04_fragments/ring_shared_atoms.csv",
     "every row of the file",
     len),
    ("body_atoms_held_in_a_ring_the_core_covers_in_part",
     "04_fragments/ring_shared_atoms.csv",
     "the sum of ring_shared_count over the file",
     lambda rows: sum(_int(r, "ring_shared_count") for r in rows)),

    ("definitions_offered",
     "05_single_edits/operator_definition_space.csv",
     "the sum of offered over the file",
     lambda rows: sum(_int(r, "offered") for r in rows)),
    ("definitions_surviving_the_valence_ledger",
     "05_single_edits/operator_definition_space.csv",
     "the sum of surviving_after_ledger over the file",
     lambda rows: sum(_int(r, "surviving_after_ledger") for r in rows)),
    ("definitions_surviving_the_bond_environment_condition",
     "05_single_edits/operator_definition_space.csv",
     "the sum of surviving_after_bond_environment over the file",
     lambda rows: sum(_int(r, "surviving_after_bond_environment")
                      for r in rows)),
    ("definitions_surviving_the_atom_environment_condition",
     "05_single_edits/operator_definition_space.csv",
     "the sum of surviving_after_atom_environment over the file",
     lambda rows: sum(_int(r, "surviving_after_atom_environment")
                      for r in rows)),
    ("definitions_accepted",
     "05_single_edits/operator_definition_space.csv",
     "the sum of accepted over the file",
     lambda rows: sum(_int(r, "accepted") for r in rows)),
    ("definitions_refused_by_the_terminal_halogen_restriction",
     "05_single_edits/operator_definition_space.csv",
     "the sum of terminal_halogen over the file",
     lambda rows: sum(_int(r, "terminal_halogen") for r in rows)),
    ("definitions_refused_for_being_above_the_lowest_valence_with_hydrogen",
     "05_single_edits/operator_definition_space.csv",
     "the sum of above_lowest_valence_with_hydrogen over the file",
     lambda rows: sum(_int(r, "above_lowest_valence_with_hydrogen")
                      for r in rows)),

    ("examinations",
     "05_single_edits/examination_records.csv",
     "every row of the file",
     len),
    ("examinations_at_the_reported_allowance",
     "05_single_edits/examination_records.csv",
     "rows whose order_index is at most 50",
     lambda rows: sum(1 for r in rows if _int(r, "order_index") <= 50)),
    ("products_certified_at_the_reported_allowance",
     "05_single_edits/examination_records.csv",
     "rows whose order_index is at most 50 and whose certified flag is 1",
     lambda rows: sum(1 for r in rows if _int(r, "order_index") <= 50
                      and r["certified"] == "1")),
    ("distinct_products_at_the_reported_allowance",
     "05_single_edits/examination_records.csv",
     "distinct product_smiles among certified rows at order_index at most 50",
     lambda rows: len({r["product_smiles"] for r in rows
                       if _int(r, "order_index") <= 50
                       and r["certified"] == "1" and r["product_smiles"]})),
    ("products_identical_to_the_parent",
     "05_single_edits/examination_records.csv",
     "rows whose identical_to_parent flag is 1",
     lambda rows: sum(1 for r in rows
                      if r.get("identical_to_parent") == "1")),
    ("products_refused_by_certification",
     "05_single_edits/refusal_records.csv",
     "every row of the file, which is expected to be empty",
     len),
    ("median_similarity_between_a_product_and_the_ion_it_came_from",
     "05_single_edits/examination_records.csv",
     "the median of similarity_to_parent over certified rows, unrounded",
     lambda rows: statistics.median(
         [v for v in (_float(r, "similarity_to_parent") for r in rows
                      if r["certified"] == "1") if v is not None] or [0.0])),
    ("judge_disagreements_with_certification",
     "05_single_edits/examination_records.csv",
     "rows whose judge_agrees_with_certification flag is 0",
     lambda rows: sum(1 for r in rows
                      if r.get("judge_agrees_with_certification") == "0")),

    ("atoms_an_edit_touched_or_built",
     "05_single_edits/valence_ledger.csv",
     "every row of the file, one per atom an edit touched or built",
     len),
    ("atoms_kept_as_they_were_and_merely_touched",
     "05_single_edits/valence_ledger.csv",
     "rows whose role is retained and touched: the edit changed what their "
     "bonds spend and nothing else about them",
     lambda rows: sum(1 for r in rows
                      if r["role"] == "retained and touched")),
    ("of_those_the_number_whose_total_valence_was_conserved",
     "05_single_edits/valence_ledger.csv",
     "retained and touched rows whose total valence before equals total "
     "valence after. The ledger is one subtraction and it requires this of "
     "every one of them, so any shortfall is a defect and not a result",
     lambda rows: sum(1 for r in rows
                      if r["role"] == "retained and touched"
                      and r["total_valence_before"]
                      == r["total_valence_after"])),

    ("signature_collisions",
     "09_validation/signature_collisions.csv",
     "every row of the file, which must be empty",
     len),
    ("perception_check_mismatches",
     "09_validation/perception_check.csv",
     "rows whose agrees flag is 0, which must be none",
     lambda rows: sum(1 for r in rows if r.get("agrees") == "0")),
    ("rings_traced_to_their_source",
     "09_validation/perception_check.csv",
     "every row of the file",
     len),

    ("string_control_attempts",
     "06_string_controls/string_control_records.csv",
     "every row of the file",
     len),
    ("string_control_products_accepted",
     "06_string_controls/string_control_records.csv",
     "rows whose accepted flag is 1",
     lambda rows: sum(1 for r in rows if r.get("accepted") == "1")),

    ("sequential_rounds_recorded",
     "07_sequential/trajectory_steps.csv",
     "every row of the file",
     len),
    ("sequential_candidates_recorded",
     "07_sequential/trajectory_candidates.csv",
     "every row of the file",
     len),
    ("sequential_roots",
     "07_sequential/root_selection.csv",
     "every row of the file, one per core",
     len),
)


def build(root: str) -> Dict[str, Any]:
    """Read the archive and write every quantity the paper may quote."""
    out: Dict[str, Any] = {}
    missing: List[str] = []
    cache: Dict[str, Optional[Rows]] = {}
    for name, relative, rule, how in TABLE:
        if relative not in cache:
            cache[relative] = read(root, relative)
        rows = cache[relative]
        if rows is None:
            missing.append(relative)
            out[name] = {"value": None, "file": relative, "rule": rule,
                         "note": "the file was not written by this run"}
            continue
        out[name] = {"value": how(rows), "file": relative, "rule": rule}
    return {"facts": out,
            "files_missing": sorted(set(missing)),
            "note": "Every value above was computed by reading the file named "
                    "beside it. The paper quotes this table and nothing else."}
