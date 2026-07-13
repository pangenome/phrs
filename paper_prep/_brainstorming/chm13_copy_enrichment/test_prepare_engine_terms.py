#!/usr/bin/env python3

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import prepare_engine_terms as bridge


def write(path, text):
    path.write_text(text)
    return path


class BridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.loci = write(self.root / "loci.tsv",
            "locus_id\tgene_name\tchromosome\tstart1\tend1\tstrand\n"
            "G1\tDUP\tchr1\t1\t10\t+\nG2\tDUP\tchr2\t21\t30\t-\n")
        self.universe = write(self.root / "universe.tsv",
            "copy_id\tseqid\tstart_1based\tend_1based_inclusive\tstrand\tgff_id\n"
            "C1\tchr1\t1\t10\t+\tG1\nC2\tchr2\t21\t30\t-\tG2\n")
        self.metadata = write(self.root / "metadata.tsv",
            "source\tterm_id\tterm_name\tnamespace\n"
            "GO\tGO:1\tprocess\tbiological_process\nbiotype\tBT:x\ttype x\tgene_biotype\n")
        self.edges = write(self.root / "edges.tsv",
            "copy_id\tsource\tterm_id\nC1\tGO\tGO:1\nC2\tGO\tGO:1\nC1\tbiotype\tBT:x\n")

    def tearDown(self):
        self.temp.cleanup()

    def test_coordinate_bridge_retains_duplicate_symbol_copies(self):
        manifest = bridge.build_collections(self.loci, self.universe, self.edges,
                                            self.metadata, self.root / "out")
        self.assertEqual(manifest["physical_loci"], 2)
        self.assertEqual(manifest["source_edges"], 3)
        with gzip.open(self.root / "out/GO_BP.tsv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["locus_id"] for row in rows}, {"G1", "G2"})
        self.assertEqual({row["copy_id"] for row in rows}, {"C1", "C2"})
        self.assertEqual({row["frozen_source"] for row in rows}, {"GO"})

    def test_coordinate_mismatch_is_a_hard_error(self):
        bad = write(self.root / "bad.tsv",
            "copy_id\tseqid\tstart_1based\tend_1based_inclusive\tstrand\tgff_id\n"
            "C1\tchr1\t2\t10\t+\tG1\nC2\tchr2\t21\t30\t-\tG2\n")
        with self.assertRaisesRegex(ValueError, "no exact analysis-locus"):
            bridge.build_crosswalk(self.loci, bad)

    def test_duplicate_source_edge_is_rejected(self):
        duplicate = write(self.root / "duplicate.tsv",
            self.edges.read_text() + "C1\tGO\tGO:1\n")
        with self.assertRaisesRegex(ValueError, "duplicate frozen"):
            bridge.build_collections(self.loci, self.universe, duplicate,
                                     self.metadata, self.root / "duplicate-out")


if __name__ == "__main__":
    unittest.main()
