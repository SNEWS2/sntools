# coding=utf-8
try:
    import __builtin__ as builtins # Python 2.7
except:
    import builtins # Python 3
import unittest

from . import es as ibd # 😅

builtins._flavor = 'e'

class CrossSectionTest(unittest.TestCase):
    # test a few sensible values
    def test_dSigma_dE(self):
        test_values = (
                       ( 3,  1.7, 2.3427673889264988e-23),
                       (20, 18.0, 2.2641725887270564e-23),
                       (20, 18.4, 2.2579850452106914e-23),
                       (50, 48.7, 2.2212735153636973e-23),
                       (100,97.0, 2.216109309506134e-23)
                      )
        for (eNu, eE, result) in test_values:
            self.assertAlmostEqual(ibd.dSigma_dE(eNu, eE), result, delta=1e-6*result)


    def test_dSigma_dE_edgecases(self):
        # test unphysical values
        self.assertEqual(ibd.dSigma_dE(ibd.bounds_eNu[0] - 1e-6, ibd.eE_min), 0) # eNu too small
        self.assertEqual(ibd.dSigma_dE(20, ibd.eE_min - 1e-6), 0) # eE too small
        self.assertEqual(ibd.dSigma_dE(20, 20.259), 0) # eE too large

    def test_dSigma_dCosT(self):
        test_values = (
                       ( 3, -0.50, 0.0),
                       ( 3,  0.00, 0.0),
                       ( 3,  0.50, 0.0),
                       (20, -0.99, 0.0),
                       (20,  0.00, 0.0),
                       (20,  0.99, 9.55532554942379e-21),
                       (50, -0.50, 0.0),
                       (50,  0.00, 0.0),
                       (50,  0.50, 4.4738155338186293e-23),
                      )
        for (eNu, cosT, result) in test_values:
            self.assertAlmostEqual(ibd.dSigma_dCosT(eNu, cosT), result, delta=1e-6*result)

    def test_dSigma_dCosT_edgecases(self):
        # unphysical values of cos(theta)
        self.assertEqual(ibd.dSigma_dCosT(20, 1.01), 0)
        self.assertEqual(ibd.dSigma_dCosT(20, -1.01), 0)

    def test_get_eE(self):
        test_values = (
                       ( 3, -0.50, 0.7391887582525531),
                       ( 3,  0.00, 0.5109989),
                       ( 3,  0.50, 0.7391887582525531),
                       (20, -0.99, 14.490377038079016),
                       (20,  0.00, 0.5109989),
                       (20,  0.99, 14.490377038079016),
                       (50, -0.50, 0.8425819323212276),
                       (50,  0.00, 0.5109989),
                       (50,  0.50, 0.8425819323212276),
                      )
        for (eNu, cosT, result) in test_values:
            self.assertAlmostEqual(ibd.get_eE(eNu, cosT), result)

    def test_bounds_eE(self):
        test_values = (
                       ( 3, 0.77, 3.275551663171255),
                       (20, 0.77, 20.258722276922214),
                       (50, 0.77, 50.25679841169901),
                      )
        for (eNu, eE_min, eE_max) in test_values:
            self.assertAlmostEqual(ibd.bounds_eE(eNu)[0], eE_min)
            self.assertAlmostEqual(ibd.bounds_eE(eNu)[1], eE_max)

    def test_bounds_eNu(self):
        self.assertAlmostEqual(ibd.bounds_eNu[0], 0.4175023400112732)
        self.assertEqual(ibd.bounds_eNu[1], 100)


if __name__ == '__main__':
    unittest.main()
