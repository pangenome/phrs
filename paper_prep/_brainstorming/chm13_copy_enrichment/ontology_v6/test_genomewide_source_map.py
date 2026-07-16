#!/usr/bin/env python3
"""Regression tests for the released V6 genome-wide source map."""

import csv
import gzip
import importlib.util
import json
import pathlib
import unittest


HERE = pathlib.Path(__file__).resolve().parent


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_module("v6_builder", HERE / "build_genomewide_source_map.py")
checker = load_module("v6_checker", HERE / "check_genomewide_source_map.py")


class HelperTests(unittest.TestCase):
    def test_named_cohort_rules_are_exact_and_alias_aware(self):
        self.assertEqual(builder.cohort_for("DUX4L6", ""), "DUX4_DUX4L")
        self.assertEqual(builder.cohort_for("DBET", "DUX4L30"), "DUX4_DUX4L")
        self.assertEqual(builder.cohort_for("TUBB8P8", ""), "TUBB8")
        self.assertEqual(builder.cohort_for("OR4F29", ""), "OR4F")
        self.assertEqual(builder.cohort_for("WASHC5-AS1", ""), "WASH")
        self.assertEqual(builder.cohort_for("NOTDUX4", ""), "")

    def test_evidence_hash_depends_on_physical_copy(self):
        row = {field: "x" for field in builder.EVIDENCE_FIELDS[:-1]}
        first = builder.row_hash(row, builder.EVIDENCE_FIELDS[:-1])
        row["copy_id"] = "y"
        second = builder.row_hash(row, builder.EVIDENCE_FIELDS[:-1])
        self.assertNotEqual(first, second)

    def test_assignment_schema_is_target_blind(self):
        self.assertFalse(any("phr" in field.lower() for field in builder.ASSIGNMENT_FIELDS))
        self.assertIn("own_annotation_id", builder.ASSIGNMENT_FIELDS)
        self.assertIn("functional_source_id", builder.ASSIGNMENT_FIELDS)
        self.assertIn("functional_mapping_status", builder.ASSIGNMENT_FIELDS)


class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = checker.check_release()
        cls.validation = json.loads(
            (HERE / "GENOMEWIDE_SOURCE_MAP_VALIDATION.json").read_text()
        )

    def test_independent_full_release_check_passes(self):
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["physical_loci"], 61312)
        self.assertEqual(self.result["phr_midpoint_copies"], 402)

    def test_every_physical_copy_and_evidence_record_is_retained(self):
        with gzip.open(str(HERE / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz"),
                       "rt", encoding="utf-8", newline="") as handle:
            assignments = list(csv.DictReader(handle, delimiter="\t"))
        with gzip.open(str(HERE / "COPY_SOURCE_EVIDENCE.tsv.gz"),
                       "rt", encoding="utf-8", newline="") as handle:
            evidence = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(assignments), 61312)
        self.assertEqual(len({row["copy_id"] for row in assignments}), 61312)
        self.assertEqual(len(evidence), 61312)
        self.assertTrue(all(row["physical_copy_cn"] == "1" for row in assignments))

    def test_named_pilot_counts_and_dux_burden(self):
        with open(str(HERE / "NAMED_COHORT_AUDIT.tsv"),
                  "rt", encoding="utf-8", newline="") as handle:
            rows = {row["cohort"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(
            (rows["DUX4_DUX4L"]["phr_physical_copies"],
             rows["DUX4_DUX4L"]["phr_ontology_contributors"]),
            ("68", "65"),
        )
        self.assertEqual(rows["DDX11L"]["phr_ontology_contributors"], "10")
        self.assertEqual(rows["TUBB8"]["phr_ontology_contributors"], "2")
        self.assertEqual(rows["OR4F"]["phr_ontology_contributors"], "4")
        self.assertEqual(rows["WASH"]["phr_ontology_contributors"], "9")
        with gzip.open(str(HERE / "NAMED_COHORT_TERM_BURDENS.tsv.gz"),
                       "rt", encoding="utf-8", newline="") as handle:
            dux = [row for row in csv.DictReader(handle, delimiter="\t")
                   if row["cohort"] == "DUX4_DUX4L"
                   and row["functional_source_symbol"] == "DUX4"]
        self.assertTrue(dux)
        self.assertTrue(all(row["phr_midpoint_physical_copy_burden"] == "65" for row in dux))

    def test_phr_join_is_recorded_after_source_freeze(self):
        self.assertTrue(self.validation["source_assignment_frozen_before_phr_open"])
        self.assertEqual(self.validation["phr_join_stage"], "after_source_assignment_freeze")
        with open(str(HERE / "INPUT_MANIFEST.tsv"),
                  "rt", encoding="utf-8", newline="") as handle:
            manifest = {row["role"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(manifest["phr_bed"]["stage"], "after_source_freeze_target_join")

    def test_release_is_mapping_only(self):
        self.assertFalse(self.validation["enrichment_run"])
        self.assertEqual(self.validation["scope"], "copy_to_ontology_mapping_only_no_enrichment")
        names = {path.name.lower() for path in HERE.iterdir()}
        self.assertFalse(any("enrichment" in name or "pvalue" in name for name in names))


if __name__ == "__main__":
    unittest.main()
