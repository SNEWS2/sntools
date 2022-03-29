import unittest

from sntools.interaction_channels import o16eb
from ._crosssectiontest import CrossSectionTest


class O16EBTest(CrossSectionTest):
    c = o16eb.Channel('eb')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (25, 13.77, 6.512218605324598e-19),
        (25, 13.7709, 6.512218605324598e-19),
        (50, 38.77, 1.2264310652357698e-17),
        (50, 31.50, 1.4897566364999385e-17),
        (50, 28.46, 9.710356625778577e-18),
        (50, 24.62, 9.551805523349172e-18),
        (75, 56.50, 6.83411669272891e-17),
        (75, 53.46, 4.659646123960687e-17),
        (100, 88.77, 1.0705158022625537e-16),
        (100, 74.62, 1.2949603633078018e-16),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (11.20, 0, 0),  # eNu too small
        (40, 14.6189, 0),  # eE too small
        (40, 28.7711, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (20, -0.99, 1.3333149289715016),
        (20, 0.00, 1.0),
        (20, 0.99, 0.6666850710284984),
        (50, -0.50, 1.3861557655727768),
        (50, 0.00, 1.0),
        (50, 0.50, 0.6138442344272234),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (18, -0.99, 6.77),
        (18, 0.00, 6.77),
        (18, 0.99, 6.77),
        # testing higher energies is more complicated, since there are multiple possible return values
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (20, 1.499, 8.771),
        (50, 24.619, 38.771),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 12.00


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
