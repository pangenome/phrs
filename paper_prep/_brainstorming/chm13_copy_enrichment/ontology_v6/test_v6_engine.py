#!/usr/bin/env python3
"""Regression tests for the V6 physical-CN ontology runner."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import v6_engine as V


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class V6EngineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_tsv(self, name: str, fields: list[str], rows: list[dict[str, object]]) -> Path:
        path = self.root / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def minimal_gate(self) -> tuple[Path, Path]:
        directory = self.root / "gate"
        repo = self.root / "repo"
        directory.mkdir()
        repo.mkdir()
        artifact = directory / "artifact.tsv"
        artifact.write_text("frozen\n", encoding="utf-8")
        audited = directory / V.AUDITED_MANIFEST.name
        self.write_manifest(audited, [{
            "path": artifact.name, "bytes": artifact.stat().st_size, "sha256": digest(artifact),
        }], ["path", "bytes", "sha256"])
        inputs = directory / V.INPUT_MANIFEST.name
        self.write_manifest(inputs, [], ["role", "stage", "path", "bytes", "sha256", "expected_sha256"])
        gate = {
            "status": "PASS", "enrichment_authorized": True, "enrichment_run": False,
            "audited_release_manifest_sha256": digest(audited),
            "physical_copy_cn": V.EXPECTED_PHYSICAL_CN,
            "ontology_eligible_copies": V.EXPECTED_ELIGIBLE_CN,
            "physical_copy_term_edges": V.EXPECTED_EDGE_CN,
            "direct_physical_copy_term_edges": V.EXPECTED_DIRECT_EDGE_CN,
            "phr_midpoint_copies": V.EXPECTED_PHR_MIDPOINT_CN,
            "phr_any_overlap_copies": V.EXPECTED_PHR_OVERLAP_CN,
            "phr_midpoint_copy_term_edges": V.EXPECTED_PHR_MIDPOINT_EDGE_CN,
            "named_cohorts": {
                name: {"phr_ontology_contributors": count}
                for name, count in V.NAMED_EXPECTED.items()
            },
        }
        (directory / V.GATE_JSON.name).write_text(json.dumps(gate), encoding="utf-8")
        return directory, repo

    @staticmethod
    def write_manifest(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_gate_passes_exact_bytes_and_refuses_digest_drift(self) -> None:
        directory, repo = self.minimal_gate()
        self.assertTrue(V.validate_gate(directory, repo)["pass"])
        (directory / "artifact.tsv").write_text("changed\n", encoding="utf-8")
        report = V.validate_gate(directory, repo)
        self.assertFalse(report["pass"])
        self.assertTrue(any(value.startswith("audited_") for value in report["errors"]))

    def test_gate_refuses_nonpass_even_when_checksums_match(self) -> None:
        directory, repo = self.minimal_gate()
        gate_path = directory / V.GATE_JSON.name
        gate = json.loads(gate_path.read_text())
        gate["status"] = "BLOCK"
        gate_path.write_text(json.dumps(gate), encoding="utf-8")
        report = V.validate_gate(directory, repo)
        self.assertFalse(report["pass"])
        self.assertIn("gate_status:BLOCK", report["errors"])

    def synthetic_freeze_inputs(self, target_column: bool = False) -> tuple[Path, Path]:
        fields = ["functional_source_id", "physical_copy_cn", "seqid", "start0", "end0"]
        if target_column:
            fields.append("phr_midpoint_cn")
        assignment_rows = [
            {"functional_source_id": "S1", "physical_copy_cn": 3, "seqid": "chr1",
             "start0": 1, "end0": 2},
            {"functional_source_id": "S2", "physical_copy_cn": 2, "seqid": "chr2",
             "start0": 11, "end0": 12},
        ]
        if target_column:
            for row in assignment_rows:
                row["phr_midpoint_cn"] = 0
        assignments = self.write_tsv("assignments.tsv", fields, assignment_rows)
        terms = self.write_tsv("terms.tsv", [
            "functional_source_id", "ontology", "relation", "namespace", "term_id", "term_name",
        ], [
            {"functional_source_id": "S1", "ontology": "GO", "relation": "direct",
             "namespace": "biological_process", "term_id": "GO:1", "term_name": "one"},
            {"functional_source_id": "S2", "ontology": "GO", "relation": "direct",
             "namespace": "biological_process", "term_id": "GO:1", "term_name": "one"},
            {"functional_source_id": "S1", "ontology": "GO", "relation": "ancestor",
             "namespace": "biological_process", "term_id": "GO:1", "term_name": "one"},
        ])
        return assignments, terms

    def test_hypothesis_freeze_sums_cn_and_separates_direct_ancestor(self) -> None:
        assignments, terms = self.synthetic_freeze_inputs()
        rows = V.build_hypothesis_rows(assignments, terms, {"chr1": (10, 20), "chr2": (10, 20)})
        self.assertEqual([(row["relation"], row["genome_physical_copy_burden"]) for row in rows],
                         [("direct", 5), ("ancestor", 3)])
        self.assertEqual(rows[0]["genome_arm_count"], 2)
        self.assertEqual(rows[0]["multiplicity_family_id"], "GO_BP")

    def test_hypothesis_freeze_refuses_target_columns(self) -> None:
        assignments, terms = self.synthetic_freeze_inputs(target_column=True)
        with self.assertRaisesRegex(ValueError, "target-membership"):
            V.build_hypothesis_rows(assignments, terms, {"chr1": (10, 20), "chr2": (10, 20)})

    def test_weighted_recount_is_sum_of_physical_cn_not_unique_labels(self) -> None:
        arrays = V.WeightedTermArrays([
            np.asarray([0, 1], dtype=np.int32),
            np.asarray([0], dtype=np.int32),
            np.empty(0, dtype=np.int32),
        ], np.asarray([65, 2, 9], dtype=np.uint32))
        observed = V.recount(np.asarray([0, 1, 2]), arrays, 2)
        np.testing.assert_array_equal(observed, [67, 65])
        self.assertNotEqual(int(observed[0]), 2)  # never a unique-row/source count
        self.assertNotEqual(int(observed[1]), 1)  # DUX-like CN is retained

    def hypothesis(self, index: int, collection: str, relation: str = "direct") -> V.Hypothesis:
        return V.Hypothesis(index, collection, relation, "GO", "biological_process",
                            f"GO:{index}", f"term {index}", 5, 2, 1)

    def test_multiplicity_families_combine_relation_layers_by_collection(self) -> None:
        rows = [
            self.hypothesis(0, "GO_BP", "direct"), self.hypothesis(1, "GO_BP", "ancestor"),
            self.hypothesis(2, "GO_MF"), self.hypothesis(3, "GO_CC"),
            V.Hypothesis(4, "Reactome", "direct", "Reactome", "pathway", "R-HSA-4", "four", 5, 2, 1),
        ]
        self.assertEqual(V.family_ranges(rows), [
            ("GO_BP", 0, 2), ("GO_MF", 2, 3), ("GO_CC", 3, 4), ("Reactome", 4, 5),
        ])

    def test_stage_screen_never_selectively_extends_maxT_unresolved_row(self) -> None:
        row = {
            "assignment": "midpoint", "global_maxT_sequential_lower": 0.04,
            "global_maxT_sequential_upper": 0.06, "collection_maxT_sequential_lower": 0.01,
            "collection_maxT_sequential_upper": 0.02, "bh_sequential_lower": 0.01,
            "bh_sequential_upper": 0.06, "raw_exceedances": 0, "p_empirical": 1e-5,
            "bh_q": 0.01, "observed_physical_copy_burden": 65, "null_q975": 2, "null_max": 3,
        }
        self.assertEqual(V.candidate_reasons(row), ["MAXT_UNRESOLVED_NO_SELECTIVE_EXTENSION"])

    def test_stage_screen_extends_only_bh_unresolved_after_maxT_resolves(self) -> None:
        row = {
            "assignment": "midpoint", "global_maxT_sequential_lower": 0.001,
            "global_maxT_sequential_upper": 0.01, "collection_maxT_sequential_lower": 0.001,
            "collection_maxT_sequential_upper": 0.01, "bh_sequential_lower": 0.001,
            "bh_sequential_upper": 0.06, "raw_exceedances": 12, "p_empirical": 0.00013,
            "bh_q": 0.01, "observed_physical_copy_burden": 65, "null_q975": 2, "null_max": 3,
        }
        reasons = V.candidate_reasons(row)
        self.assertIn("BH_INTERVAL_STRADDLES_0.05", reasons)
        self.assertIn("RAW_EXCEEDANCES_BELOW_100", reasons)

    def test_stage_screen_does_not_extend_resolved_low_tail(self) -> None:
        row = {
            "assignment": "midpoint", "global_maxT_sequential_lower": 0.001,
            "global_maxT_sequential_upper": 0.01, "collection_maxT_sequential_lower": 0.001,
            "collection_maxT_sequential_upper": 0.01, "bh_sequential_lower": 0.001,
            "bh_sequential_upper": 0.02, "raw_exceedances": 2, "p_empirical": 0.00003,
            "bh_q": 0.01, "observed_physical_copy_burden": 65, "null_q975": 2, "null_max": 3,
        }
        self.assertEqual(V.candidate_reasons(row), [])

    def test_stage_screen_does_not_call_maxT_unresolved_when_bh_is_nonpass(self) -> None:
        row = {
            "assignment": "midpoint", "global_maxT_sequential_lower": 0.04,
            "global_maxT_sequential_upper": 0.06, "collection_maxT_sequential_lower": 0.04,
            "collection_maxT_sequential_upper": 0.06, "bh_sequential_lower": 0.2,
            "bh_sequential_upper": 0.3, "raw_exceedances": 50_000, "p_empirical": 0.5,
            "bh_q": 0.5, "observed_physical_copy_burden": 1, "null_q975": 2,
            "null_max": 3,
        }
        self.assertEqual(V.classify_primary(row), "CERTIFIED_NONPASS")
        self.assertEqual(V.candidate_reasons(row), [])

    def test_real_sampler_preserves_all_rigid_block_invariants(self) -> None:
        arms, _intervals, _loci, _genome, blocks = V.N.engine_objects()
        sampler = V.E.RegionSampler(blocks, arms, V.E.load_masks(None, arms), "primary", min_candidates=100)
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(37)))
        for replicate in range(1, 11):
            placements = sampler.sample(rng, replicate)
            self.assertEqual(len(placements), 37)
            for placed in placements:
                source = next(row for row in blocks if row.block_id == placed.block_id)
                self.assertEqual(placed.arm, source.source_arm)
                self.assertEqual(placed.components, source.components)
                self.assertEqual(placed.end - placed.start, source.span)
                midpoint = placed.start + source.midpoint_offset
                self.assertEqual(V.E.stratum_index(V.E.terminal_distance(arms[placed.arm], midpoint)),
                                 source.stratum)


if __name__ == "__main__":
    unittest.main()
