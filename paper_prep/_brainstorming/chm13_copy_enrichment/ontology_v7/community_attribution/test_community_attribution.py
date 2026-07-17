#!/usr/bin/env python3
"""Tests for the V7 post-inference community-attribution release."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("build_community_attribution", HERE / "build_community_attribution.py")
assert SPEC and SPEC.loader
attribution = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(attribution)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


class CommunityAttributionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temp = tempfile.TemporaryDirectory()
        cls.generated = Path(cls.temp.name)
        cls.validation = attribution.generate(cls.generated)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_release_validation_passes_all_gates(self) -> None:
        self.assertEqual(self.validation["status"], "PASS", self.validation["errors"])
        self.assertEqual(self.validation["checks_passed"], self.validation["checks_total"])
        self.assertEqual(self.validation["primary_supported_exact_terms"], 209)
        self.assertEqual(self.validation["overlap_sensitivity_supported_exact_terms"], 211)

    def test_exact_rows_preserve_copy_ids_and_derive_dominance(self) -> None:
        rows = read_tsv(self.generated / "COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv")
        self.assertTrue(rows)
        for row in rows:
            copy_ids = json.loads(row["copy_ids_json"])
            self.assertEqual(len(copy_ids), int(row["coordinate_copy_count"]))
            self.assertEqual(len(copy_ids), int(row["community_contributor_copy_burden"]))
            self.assertEqual(int(row["dominant_community_flag"]),
                             int(float(row["community_fraction_of_a_T"]) >= 0.5))
            self.assertTrue(all("|" in copy_id for copy_id in copy_ids))
        headline = [row for row in rows if row["analysis_layer"] == "HEADLINE_PRIMARY_MIDPOINT"]
        sensitivity = [row for row in rows if row["analysis_layer"] == "SENSITIVITY_ANY_OVERLAP"]
        self.assertTrue(all(row["support_status"] == "PRIMARY_SUPPORTED" and row["assignment"] == "midpoint"
                            for row in headline))
        self.assertTrue(all(row["support_status"] == "SENSITIVITY_SUPPORTED" and row["assignment"] == "overlap"
                            for row in sensitivity))

    def test_functional_classes_have_no_inferential_columns(self) -> None:
        prohibited = {"p_value", "p_exact_upper", "bh_q_within_collection", "by_q_within_collection",
                      "holm_p_global", "bonferroni_p_global", "primary_support", "sensitivity_support"}
        for name in ("FUNCTIONAL_CLASS_DEFINITIONS.tsv", "COMMUNITY_FUNCTIONAL_CLASSES.tsv",
                     "FUNCTIONAL_CLASS_COMMUNITY_SUMMARY.tsv", "ARM_FUNCTIONAL_CLASS_PROJECTION.tsv"):
            rows = read_tsv(self.generated / name)
            self.assertTrue(rows)
            self.assertFalse(prohibited & set(rows[0]), name)
        definitions = read_tsv(self.generated / "FUNCTIONAL_CLASS_DEFINITIONS.tsv")
        self.assertEqual(len(definitions), 6)
        self.assertTrue(all(row["definition_status"] == "PRIMARY_SUPPORTED" for row in definitions))
        community_rows = read_tsv(self.generated / "COMMUNITY_FUNCTIONAL_CLASSES.tsv")
        self.assertEqual({row["community"] for row in community_rows}, {f"C{index}" for index in range(1, 16)})
        zero_burden = {row["community"] for row in community_rows
                       if row["row_type"] == "COMMUNITY_TOTAL_NO_DEFINED_CLASS" and
                       row["community_ontology_eligible_phr_copy_burden"] == "0"}
        self.assertEqual(zero_burden, {"C4", "C10", "C13"})

    def test_arm_projection_contains_all_chromosome_arms(self) -> None:
        rows = read_tsv(self.generated / "ARM_FUNCTIONAL_CLASS_PROJECTION.tsv")
        self.assertGreater(len(rows), 48)
        labels = {row["chrom_arm_label"] for row in rows}
        expected = {f"{chrom}{arm}" for chrom in list(map(str, range(1, 23))) + ["X", "Y"] for arm in ("p", "q")}
        self.assertEqual(labels, expected)
        no_signal = {row["chrom_arm_label"] for row in rows if row["row_type"] == "NO_SIGNAL_ARM"}
        self.assertEqual(no_signal, {"2p", "3p", "5p", "8q", "11q", "14q", "18q"})
        c1_arms = {(row["chrom_arm_label"], row["community"], row["display_label"], row["class_unique_coordinate_copy_burden"])
                   for row in rows if row["class_id"] == "DUX4_ZGA_TRANSCRIPTION_NUCLEAR_ENVELOPE_CELL_CYCLE"}
        self.assertEqual(c1_arms, {
            ("4q", "C1", "DUX4 / ZGA / transcription / nuclear envelope / cell cycle", "33"),
            ("10q", "C1", "DUX4 / ZGA / transcription / nuclear envelope / cell cycle", "32"),
        })
        report = (self.generated / "COMMUNITY_FUNCTIONAL_ATTRIBUTION_REPORT.md").read_text(encoding="utf-8")
        self.assertIn("## Chromosome-arm flashcards", report)
        self.assertIn("**4q** — C1", report)
        self.assertIn("**2p** — NO_SIGNAL_ARM", report)

    def test_committed_release_is_reproducible(self) -> None:
        names = (
            "COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv",
            "COMMUNITY_SOURCE_COPY_SUMMARY.tsv",
            "FUNCTIONAL_CLASS_DEFINITIONS.tsv",
            "COMMUNITY_FUNCTIONAL_CLASSES.tsv",
            "FUNCTIONAL_CLASS_COMMUNITY_SUMMARY.tsv",
            "ARM_FUNCTIONAL_CLASS_PROJECTION.tsv",
            "COMMUNITY_FUNCTIONAL_ATTRIBUTION_REPORT.md",
            "COMMUNITY_FUNCTIONAL_ATTRIBUTION_VALIDATION.json",
            "OUTPUT_MANIFEST.sha256.tsv",
        )
        for name in names:
            self.assertEqual((HERE / name).read_bytes(), (self.generated / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
