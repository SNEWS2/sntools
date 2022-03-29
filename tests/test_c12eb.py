import unittest

from sntools.interaction_channels import c12eb
from ._crosssectiontest import CrossSectionTest


class C12EBTest(CrossSectionTest):
    c = c12eb.Channel('eb')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (25, 10.61, 2.63032746305625e-18),
        (25, 10.6109, 2.63032746305625e-18),
        (50, 35.61, 5.493977104481186e-17),
        (75, 60.61, 2.0634285413474464e-16),
        (100, 85.61, 4.993427358604262e-16),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (14.389, 0, 0),  # eNu too small
        (25, 10.6089, 0),  # eE too small
        (25, 10.6111, 0),  # eE too large
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
        (20, -0.99, 5.61),
        (20, 0.00, 5.61),
        (20, 0.99, 5.61),
        # testing higher energies is more complicated, since there are multiple possible return values
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (20, 5.609, 5.611),
        (50, 35.609, 35.611),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 15.16


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
