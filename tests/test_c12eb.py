import unittest

from sntools.interaction_channels import c12eb
from ._crosssectiontest import CrossSectionTest


class C12EBTest(CrossSectionTest):
    c = c12eb.Channel('eb')  # ensure we can access interaction channel module as self.c
    
    # To convert cm^2 to MeV^-2: http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    cm2_to_MeVm2 = (5.067731E10)**2  # convert cm^2 to MeV^-2: http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    e_thr = 14.39  # energy threshold of this reaction (arXiv:1507.05613, p. 74)
    epsilon = 0.001  # Used for approximating DiracDelta distribution.
    deriv_factor = 1.0/(2.0 * epsilon)
    #Factor to convert from sigma in cm^2 to dSigma/dE in MeV^-3.
    conv_factor = deriv_factor * cm2_to_MeVm2

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    #Test values from Table V in arXiv:nucl-th/9903022
    test_dSigma_dE_values = (
        (20, 20-e_thr, 0.8e-42 * conv_factor),
        (20, 20-e_thr+0.9*epsilon, 0.8e-42 * conv_factor),
        (50, 50-e_thr, 41.4e-42 * conv_factor),
        (70, 70-e_thr, 130e-42 * conv_factor),
        (100, 100-e_thr, 386e-42 * conv_factor),
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
