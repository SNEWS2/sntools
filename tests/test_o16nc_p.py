import unittest

from sntools.interaction_channels import o16nc_p
from ._crosssectiontest import CrossSectionTest


class O16ETest(CrossSectionTest):
    c = o16nc_p.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (20, 0, 3.955012213207594e-20),
        (20, 10, 3.955012213207594e-20),
        (30, 10, 9.96457622548407e-19),
        (50, 10, 1.951824209115436e-17),
        (75, 10, 1.0277407236510741e-16),
        (100, 10, 2.2985298252083094e-16),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (13.99, 0, 0),  # eNu too small
        # dSigma_dE(eNu, eE) is independent of eE, so no edge cases for eE exist
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (20, -0.99, 3.955012213207594e-23),
        (20, 0, 3.955012213207594e-23),
        (20, 0.99, 3.955012213207594e-23),
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
    test_bounds_eNu_minvalue = 14.0


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
