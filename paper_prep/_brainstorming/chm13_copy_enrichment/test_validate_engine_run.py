#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

import numpy as np

import COPY_engine
import validate_engine_run


def write(path, text):
    path.write_text(text)
    return path


class IndependentRunValidatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.arms = write(self.root / "arms.tsv",
            "arm\tchromosome\tstart0\tend0\nchr1_p\tchr1\t0\t20\n")
        self.intervals = write(self.root / "intervals.tsv",
            "phr_id\tchromosome\tarm\tstart0\tend0\nP1\tchr1\tchr1_p\t0\t4\n")
        self.loci = write(self.root / "loci.tsv",
            "locus_id\tgene_name\tchromosome\tarm\tstart0\tend0\tmidpoint0\tgene_biotype\n"
            "L1\tDUP\tchr1\tchr1_p\t0\t2\t1\tprotein_coding\n"
            "L2\tDUP\tchr1\tchr1_p\t2\t4\t3\tprotein_coding\n"
            "L3\tX\tchr1\tchr1_p\t10\t13\t11\tlncRNA\n")
        self.terms = write(self.root / "terms.tsv",
            "locus_id\tterm_id\tterm_name\tfrozen_source\tcopy_id\n"
            "L1\tT1\tterm one\tGO\tCOPY1\nL2\tT1\tterm one\tGO\tCOPY2\n"
            "L3\tT2\tterm two\tGO\tCOPY3\n")
        self.run = self.root / "run"
        arguments = ["--arms", str(self.arms), "--intervals", str(self.intervals),
                     "--loci", str(self.loci), "--terms", "GO_BP=%s" % self.terms,
                     "--min-term-loci", "1", "--min-term-arms", "1",
                     "--min-candidates", "1", "--permutations", "11", "--batch-size", "4",
                     "--seed", "72", "--allow-pilot", "--output", str(self.run)]
        self.assertEqual(COPY_engine.main(arguments), 0)

    def tearDown(self):
        self.temp.cleanup()

    def test_black_box_recount_and_inference_agree(self):
        report = validate_engine_run.validate(self.run)
        self.assertEqual(report["physical_loci"], 3)
        self.assertEqual(report["result_rows"], 12)
        self.assertTrue(all(value == "pass" for key, value in report["checks"].items()
                            if key != "impossible_contingencies"))

    def test_corrupt_cached_count_is_detected(self):
        path = next((self.run / "batches").glob("midpoint.GO_BP.counts.*.npy"))
        value = np.load(path)
        value[0, 0] += 1
        with path.open("wb") as handle:
            np.save(handle, value)
        with self.assertRaises(AssertionError):
            validate_engine_run.validate(self.run)


if __name__ == "__main__":
    unittest.main()
