import unittest

from sntools.interaction_channels import c12e
from ._crosssectiontest import CrossSectionTest


class C12ETest(CrossSectionTest):
    c = c12e.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (25, 7.662, 3.3139319366555277e-18),
        (25, 7.6629, 3.3139319366555277e-18),
        (50, 32.662, 7.669220668928699e-17),
        (75, 57.662, 4.196362653200819e-16),
        (100, 82.662, 1.2581902241845584e-15),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (17.337, 0, 0),  # eNu too small
        (25, 7.6609, 0),  # eE too small
        (25, 7.6631, 0),  # eE too large
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
        (20, -0.99, 2.662),
        (20, 0.00, 2.662),
        (20, 0.99, 2.662),
        # testing higher energies is more complicated, since there are multiple possible return values
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (20, 2.661, 2.663),
        (50, 32.661, 32.663),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 18.108


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
