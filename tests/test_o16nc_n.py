import unittest

from sntools.interaction_channels import o16nc_n
from ._crosssectiontest import CrossSectionTest


class O16ETest(CrossSectionTest):
    c = o16nc_n.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (20, 0, 2.478303107626836e-22),
        (20, 10, 2.478303107626836e-22),
        (30, 10, 1.6564823879992844e-19),
        (50, 10, 5.7013812424161415e-18),
        (75, 10, 3.37836296864102e-17),
        (100, 10, 7.845819682694285e-17),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (19.49, 0, 0),  # eNu too small
        # dSigma_dE(eNu, eE) is independent of eE, so no edge cases for eE exist
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (20, -0.99, 2.478303107626836e-25),
        (20, 0, 2.478303107626836e-25),
        (20, 0.99, 2.478303107626836e-25),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        # cannot test this, since return value is uncertain
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        # cannot test this, since return value is uncertain
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 19.5


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
