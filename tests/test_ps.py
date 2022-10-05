# coding=utf-8
import unittest

from sntools.interaction_channels import ep
from ._crosssectiontest import CrossSectionTest


class EPTest(CrossSectionTest):
    c = ep.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (3, 938.29, 3.9490477504e-20),
        (20, 938.5, 2.5754820771e-20),
        (20, 939, 3.72685198102e-20),
        (50, 943, 3.69517005422e-20),
        (100, 955, 3.37000485374e-20),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (20, c.eE_min - 1e-6, 0),  # eE too small
        (20, 20.259, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (3, -0.50, 0.0),
        (3, 0.00, 0.0),
        (3, 0.50, 4.90727709073e-22),
        (20, -0.99, 0.0),
        (20, 0.00, 0.0),
        (20, 0.99, 0.0), 
        (50, -0.50, 0.0),
        (50, 0.00, 0.0),
        (50, 0.50, 1.33894481498e-19),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (3, -0.50, 938.27684205),
        (3, 0.00, 938.272046),
        (3, 0.50, 938.27684205),
        (20, -0.99, 939.107709818),
        (20, 0.00, 938.272046),
        (20, 0.99, 939.107709818),
        (50, -0.50, 939.604282216),
        (50, 0.00, 938.272046),
        (50, 0.50, 939.604282216),
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (3, 938.272046, 938.291108303),
        (20, 938.272046, 939.089814435),
        (50, 938.272046, 943.087739555),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 0.0


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
