#!/usr/bin/env python3
"""Regression tests for the independent V6 pre-inference gate."""

import importlib.util
import pathlib
import unittest


HERE = pathlib.Path(__file__).resolve().parent


def load_gate():
    path = HERE / "independent_gate_genomewide_source_map.py"
    spec = importlib.util.spec_from_file_location("independent_v6_gate", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = load_gate()


class SourceEvidencePolicyTests(unittest.TestCase):
    def test_exact_self_requires_frozen_terms(self):
        allowed = gate.adjudicate_source(
            own_gene_id="1",
            source_gene_id="1",
            disposition="EXACT_SELF",
            relationship="exact registry identity",
            ambiguity="none",
            term_bearing_gene_ids={"1"},
            directed_relations=set(),
        )
        termless = gate.adjudicate_source(
            own_gene_id="1",
            source_gene_id="1",
            disposition="EXACT_SELF",
            relationship="exact registry identity",
            ambiguity="none",
            term_bearing_gene_ids=set(),
            directed_relations=set(),
        )
        self.assertEqual(allowed, (True, "exact_term_bearing_self"))
        self.assertEqual(termless, (False, "termless_self"))

    def test_name_family_alias_and_one_way_sequence_cannot_broadcast(self):
        for label in ("name", "family", "alias", "one_way_sequence"):
            with self.subTest(label=label):
                observed = gate.adjudicate_source(
                    own_gene_id="10",
                    source_gene_id="20",
                    disposition=label,
                    relationship=label,
                    ambiguity="none",
                    term_bearing_gene_ids={"20"},
                    directed_relations=set(),
                )
                self.assertEqual(observed, (False, "missing_exact_directed_relation"))

    def test_reverse_only_relation_is_rejected(self):
        observed = gate.adjudicate_source(
            own_gene_id="10",
            source_gene_id="20",
            disposition="EXPLICIT_RELATED_FUNCTIONAL_GENE",
            relationship="Related functional gene",
            ambiguity="none",
            term_bearing_gene_ids={"20"},
            directed_relations={("20", "10")},
        )
        self.assertEqual(observed, (False, "missing_exact_directed_relation"))

    def test_ambiguity_fails_closed_even_with_directed_relation(self):
        observed = gate.adjudicate_source(
            own_gene_id="10",
            source_gene_id="20",
            disposition="EXPLICIT_RELATED_FUNCTIONAL_GENE",
            relationship="Related functional gene",
            ambiguity="MULTIPLE_ONTOLOGY_BEARING_RELATED_FUNCTIONAL_GENES",
            term_bearing_gene_ids={"20"},
            directed_relations={("10", "20")},
        )
        self.assertEqual(observed, (False, "ambiguous"))

    def test_exact_directed_relation_is_admissible(self):
        observed = gate.adjudicate_source(
            own_gene_id="10",
            source_gene_id="20",
            disposition="EXPLICIT_RELATED_FUNCTIONAL_GENE",
            relationship="Related functional gene",
            ambiguity="none",
            term_bearing_gene_ids={"20"},
            directed_relations={("10", "20")},
        )
        self.assertEqual(observed, (True, "exact_directed_related_functional_gene"))


class CopyNumberTests(unittest.TestCase):
    def test_repeated_source_retains_every_physical_contributor(self):
        copies = [
            {"copy_id": "copy1", "source": "S", "physical": 1, "phr": 1},
            {"copy_id": "copy2", "source": "S", "physical": 1, "phr": 1},
            {"copy_id": "copy3", "source": "S", "physical": 1, "phr": 0},
            {"copy_id": "copy4", "source": "Z", "physical": 1, "phr": 0},
        ]
        burdens = gate.copy_weighted_burdens(
            copies, {"S": ("T1", "T2"), "Z": ("T3",)}
        )
        self.assertEqual(burdens[("S", "T1")], (3, 2, 3, 2))
        self.assertEqual(burdens[("S", "T2")], (3, 2, 3, 2))
        self.assertEqual(burdens[("Z", "T3")], (1, 0, 1, 0))


class ReleasedGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = gate.run_gate(write_outputs=False)

    def test_independent_release_reconstruction_passes(self):
        self.assertEqual(self.result["status"], "PASS")
        self.assertFalse(self.result["production_builder_imported"])
        self.assertFalse(self.result["enrichment_run"])
        self.assertTrue(self.result["enrichment_authorized"])

    def test_raw_copy_and_phr_universes_are_complete(self):
        self.assertEqual(self.result["physical_loci"], 61312)
        self.assertEqual(self.result["physical_copy_cn"], 61312)
        self.assertEqual(self.result["phr_midpoint_copies"], 402)
        self.assertEqual(self.result["phr_any_overlap_copies"], 412)

    def test_every_source_term_edge_and_burden_is_reconstructed(self):
        self.assertEqual(self.result["source_terms"], 1686727)
        self.assertEqual(self.result["physical_copy_term_edges"], 2929709)
        self.assertEqual(self.result["term_burden_rows"], 1686727)

    def test_reviewed_named_cohort_burdens_are_copy_counts(self):
        observed = {
            name: values["phr_ontology_contributors"]
            for name, values in self.result["named_cohorts"].items()
        }
        self.assertEqual(
            observed,
            {"DUX4_DUX4L": 65, "DDX11L": 10, "TUBB8": 2, "OR4F": 4, "WASH": 9},
        )

    def test_all_machine_checks_pass(self):
        self.assertTrue(self.result["checks"])
        self.assertTrue(all(row["status"] == "PASS" for row in self.result["checks"]))


if __name__ == "__main__":
    unittest.main()
