#!/usr/bin/env python3

import json
import unittest

import calibration_suite


class CalibrationUnitTests(unittest.TestCase):
    def test_chunk_seed_is_reproducible(self):
        first = calibration_suite.run(0, 1, 4, 39, 100, 123456)
        second = calibration_suite.run(0, 1, 4, 39, 100, 123456)
        first.pop("elapsed_seconds")
        second.pop("elapsed_seconds")
        self.assertEqual(json.dumps(first, sort_keys=True), json.dumps(second, sort_keys=True))

    def test_weighted_failure_control_detects_cluster_mismatch(self):
        import numpy as np
        rng = np.random.Generator(np.random.PCG64DXSM(991))
        result = calibration_suite.weighted_failure_control(rng, 2000)
        self.assertGreater(result["rejections"] / result["trials"], 0.15)


if __name__ == "__main__":
    unittest.main()
