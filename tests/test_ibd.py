import unittest

from sntools.interaction_channels import ibd
from ._crosssectiontest import CrossSectionTest


class IBDTest(CrossSectionTest):
    c = ibd.Channel('eb')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (3, 1.7, 6.597590600586344e-20),
        (20, 18.0, 9.756791751947533e-20),
        (20, 18.4, 9.686960987505651e-20),
        (50, 48.7, 1.003350542600863e-19),
        (100, 97.0, 9.302012628854371e-20),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (ibd.eThr - 1e-6, ibd.mE, 0),  # eNu too small
        (20, 17.841, 0),  # eE too small
        (20, 18.706, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (3, -0.50, 3.520263968746167e-22),
        (3, 0.00, 3.381936718733454e-22),
        (3, 0.50, 3.242202950723807e-22),
        (20, -0.99, 3.583761750083326e-20),
        (20, 0.00, 3.705455845604023e-20),
        (20, 0.99, 3.852469055868431e-20),
        (50, -0.50, 1.736535165001988e-19),
        (50, 0.00, 1.970150743424877e-19),
        (50, 0.50, 2.255441012065604e-19),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (3, -0.50, 1.6979103666324273),
        (3, 0.00, 1.700490693779404),
        (3, 0.50, 1.7030796578148266),
        (20, -0.99, 17.944888098364295),
        (20, 0.00, 18.31551834516354),
        (20, 0.99, 18.701794009695444),
        (50, -0.50, 45.100902584113804),
        (50, 0.00, 46.24173207227394),
        (50, 0.50, 47.44178162206098),
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (3, 1.695338635167926, 1.705677300207871),
        (20, 17.941220926704954, 18.70577898514169),
        (50, 44.01501537362355, 48.70578493879961),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 1.8060455572837861


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
