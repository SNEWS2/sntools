import unittest

from sntools.interaction_channels import c12nc
from ._crosssectiontest import CrossSectionTest


class C12NCTest(CrossSectionTest):
    c = c12nc.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (25, 15.11, 1.895519358671657e-18),
        (25, 15.1109, 1.895519358671657e-18),
        (50, 15.11, 2.3590523576277865e-17),
        (75, 15.11, 6.950951642452967e-17),
        (100, 15.11, 1.396524979034271e-16),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (15.1, 15.11, 0),  # eNu too small
        (25, 15.1089, 0),  # eE too small
        (25, 15.1111, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (20, -0.99, 0.5),
        (20, 0.00, 0.5),
        (20, 0.99, 0.5),
        (50, -0.50, 0.5),
        (50, 0.00, 0.5),
        (50, 0.50, 0.5),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (20, -0.99, 15.11),
        (20, 0.00, 15.11),
        (20, 0.99, 15.11),
        # testing higher energies is more complicated, since there are multiple possible return values
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (20, 15.109, 15.111),
        (50, 15.109, 15.111),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 15.11


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
