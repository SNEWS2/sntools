"""
CrossSectionTest

Base class for implementing tests for interaction channels.
Separate implementations for each channel must define numerical test
values as described in the comments above each test below.
"""

import unittest

# TODO: add test for `generate_event(eNu, direction)`


class CrossSectionTest(unittest.TestCase):
    # In child class, define `c` as the imported module of that channel's
    # cross-section so we can use `self.c` here.

    # In child class, define `test_dSigma_dE_values` with cross-sections
    # for some typical values of eNu, eE to run this test.
    def test_dSigma_dE(self):
        for (eNu, eE, result) in self.test_dSigma_dE_values:
            self.assertAlmostEqual(self.c.dSigma_dE(eNu, eE), result, delta=1e-6 * result)

    # In child class, define `test_dSigma_dE_edgecases_values` with
    # cross-sections for some typical values of eNu, eE to run this test.
    def test_dSigma_dE_edgecases(self):
        for (eNu, eE, result) in self.test_dSigma_dE_edgecases_values:
            self.assertEqual(self.c.dSigma_dE(eNu, eE), result)  # eE too small

    # In child class, define `test_dSigma_dCosT_values` with
    # cross-sections for some typical values of eNu, cosT to run this test.
    def test_dSigma_dCosT(self):
        for (eNu, cosT, result) in self.test_dSigma_dCosT_values:
            self.assertAlmostEqual(self.c.dSigma_dCosT(eNu, cosT), result, delta=1e-6 * result)

    # Untypical values for cosT are independent of interaction channel.
    # This test runs automatically without defining values in child class.
    def test_dSigma_dCosT_edgecases(self):
        self.assertEqual(self.c.dSigma_dCosT(20, 1.01), 0)
        self.assertEqual(self.c.dSigma_dCosT(20, -1.01), 0)

    # In child class, define `test_get_eE_values` with
    # eE for some typical values of eNu, cosT to run this test.
    def test_get_eE(self):
        for (eNu, cosT, result) in self.test_get_eE_values:
            self.assertAlmostEqual(self.c.get_eE(eNu, cosT), result)

    # In child class, define `test_bounds_eE_values` with min/max
    # values of eE for some typical values of eNu to run this test.
    def test_bounds_eE(self):
        for (eNu, eE_min, eE_max) in self.test_bounds_eE_values:
            self.assertAlmostEqual(self.c.bounds_eE(eNu)[0], eE_min)
            self.assertAlmostEqual(self.c.bounds_eE(eNu)[1], eE_max)

    # In child class, define `test_bounds_eNu_minvalue` as the minimum
    # value of eNu to run this test.
    def test_bounds_eNu(self):
        self.assertAlmostEqual(self.c.bounds_eNu[0], self.test_bounds_eNu_minvalue)
        self.assertEqual(self.c.bounds_eNu[1], 100)
