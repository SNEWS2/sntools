import unittest

from sntools.interaction_channels import o16e
from ._crosssectiontest import CrossSectionTest


class O16ETest(CrossSectionTest):
    c = o16e.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (25, 9.79, 3.851148525797014e-19),
        (25, 9.7909, 3.851148525797014e-19),
        (50, 34.79, 2.081988997028917e-17),
        (50, 27.53, 3.337300789732481e-17),
        (50, 24.49, 4.0646617371610224e-18),
        (50, 20.65, 1.809041245560864e-17),
        (75, 52.53, 2.0549061427325933e-16),
        (75, 49.49, 4.5100346438456647e-17),
        (100, 84.79, 2.9736990956760937e-16),
        (100, 70.65, 4.3495907530083017e-16),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (15.20, 0, 0),  # eNu too small
        (40, 10.6489, 0),  # eE too small
        (40, 24.7911, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (20, -0.99, 1.3303518123667377),
        (20, 0.00, 1.0),
        (20, 0.99, 0.6696481876332623),
        (50, -0.50, 1.3538353601496727),
        (50, 0.00, 1.0),
        (50, 0.50, 0.6461646398503273),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (20, -0.99, 4.79),
        (20, 0.00, 4.79),
        (20, 0.99, 4.79),
        # testing higher energies is more complicated, since there are multiple possible return values
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (20, 4.789, 4.791),
        (50, 20.649, 34.791),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 15.98


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
