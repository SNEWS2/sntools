import unittest

from . import ibd

class CrossSectionTest(unittest.TestCase):
    # test a few sensible values
    def test_dSigma_dE(self):
        test_values = (
                       ( 3,  1.7, 6.597590600586344e-20),
                       (20, 18.0, 9.756791751947533e-20),
                       (20, 18.4, 9.686960987505651e-20),
                       (50, 48.7, 1.003350542600863e-19),
                       (100,97.0, 9.302012628854371e-20)
                      )
        for (eNu, eE, result) in test_values:
            self.assertAlmostEqual(ibd.dSigma_dE(eNu, eE), result, delta=1e-6*result)


    def test_dSigma_dE_edgecases(self):
        # test unphysical values
        self.assertEqual(ibd.dSigma_dE(ibd.eThr - 1e-6, ibd.mE), 0) # eNu too small
        self.assertEqual(ibd.dSigma_dE(20, 17.941), 0) # eE too small
        self.assertEqual(ibd.dSigma_dE(20, 18.706), 0) # eE too large

    def test_dSigma_dCosT(self):
        test_values = (
                       ( 3, -0.50, 3.520263968746167e-22),
                       ( 3,  0.00, 3.381936718733454e-22),
                       ( 3,  0.50, 3.242202950723807e-22),
                       (20, -0.99, 3.583761750083326e-20),
                       (20,  0.00, 3.705455845604023e-20),
                       (20,  0.99, 3.852469055868431e-20),
                       (50, -0.50, 1.736535165001988e-19),
                       (50,  0.00, 1.970150743424877e-19),
                       (50,  0.50, 2.255441012065604e-19),
                      )
        for (eNu, cosT, result) in test_values:
            self.assertAlmostEqual(ibd.dSigma_dCosT(eNu, cosT), result, delta=1e-6*result)

    def test_dSigma_dCosT_edgecases(self):
        # unphysical values of cos(theta)
        self.assertEqual(ibd.dSigma_dCosT(20, 1.01), 0)
        self.assertEqual(ibd.dSigma_dCosT(20, -1.01), 0)

    def test_get_eE(self):
        test_values = (
                       ( 3, -0.50, 1.6979103666324273),
                       ( 3,  0.00, 1.700490693779404),
                       ( 3,  0.50, 1.7030796578148266),
                       (20, -0.99, 17.944888098364295),
                       (20,  0.00, 18.31551834516354),
                       (20,  0.99, 18.701794009695444),
                       (50, -0.50, 45.100902584113804),
                       (50,  0.00, 46.24173207227394),
                       (50,  0.50, 47.44178162206098),
                      )
        for (eNu, cosT, result) in test_values:
            self.assertAlmostEqual(ibd.get_eE(eNu, cosT), result)

    def test_bounds_eE(self):
        test_values = (
                       ( 3, 1.695338635167926, 1.705677300207871),
                       (20, 17.941220926704954, 18.70577898514169),
                       (50, 44.01501537362355, 48.70578493879961),
                      )
        for (eNu, eE_min, eE_max) in test_values:
            self.assertAlmostEqual(ibd.bounds_eE(eNu)[0], eE_min)
            self.assertAlmostEqual(ibd.bounds_eE(eNu)[1], eE_max)

    def test_bounds_eNu(self):
        self.assertAlmostEqual(ibd.eThr, 1.8060455572837861)
        self.assertEqual(ibd.bounds_eNu[0], ibd.eThr)
        self.assertEqual(ibd.bounds_eNu[1], 100)


if __name__ == '__main__':
    unittest.main()
