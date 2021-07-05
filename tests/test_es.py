# coding=utf-8
import unittest

from sntools.interaction_channels import es
from ._crosssectiontest import CrossSectionTest


class ESTest(CrossSectionTest):
    c = es.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (3, 1.7, 2.3427673889264988e-23),
        (20, 18.0, 2.2641725887270564e-23),
        (20, 18.4, 2.2579850452106914e-23),
        (50, 48.7, 2.2212735153636973e-23),
        (100, 97.0, 2.216109309506134e-23),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (c.bounds_eNu[0] - 1e-6, c.eE_min, 0),  # eNu too small
        (20, c.eE_min - 1e-6, 0),  # eE too small
        (20, 20.259, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (3, -0.50, 0.0),
        (3, 0.00, 0.0),
        (3, 0.50, 0.0),
        (20, -0.99, 0.0),
        (20, 0.00, 0.0),
        (20, 0.99, 9.55532554942379e-21),
        (50, -0.50, 0.0),
        (50, 0.00, 0.0),
        (50, 0.50, 4.4738155338186293e-23),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (3, -0.50, 0.7391887582525531),
        (3, 0.00, 0.5109989),
        (3, 0.50, 0.7391887582525531),
        (20, -0.99, 14.490377038079016),
        (20, 0.00, 0.5109989),
        (20, 0.99, 14.490377038079016),
        (50, -0.50, 0.8425819323212276),
        (50, 0.00, 0.5109989),
        (50, 0.50, 0.8425819323212276),
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (3, 0.77, 3.275551663171255),
        (20, 0.77, 20.258722276922214),
        (50, 0.77, 50.25679841169901),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 0.4175023400112732


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
