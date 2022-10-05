# coding=utf-8
import unittest

from sntools.interaction_channels import ps
from ._crosssectiontest import CrossSectionTest


class PSTest(CrossSectionTest):
    c = ps.Channel('e')  # ensure we can access interaction channel module as self.c

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_values = (
        (3, 938.29, 3.94862488664774e-20),
        (20, 938.5, 2.5754728840458534e-20),
        (20, 939, 3.7268427887733035e-20),
        (50, 943, 3.6951686796684256e-20),
        (100, 955, 3.370004553363625e-20),
    )

    # iterable with tuples (eNu, eE, dSigma_dE(eNu, eE))
    test_dSigma_dE_edgecases_values = (
        (20, c.eE_min - 1e-6, 0),  # eE too small
        (20, 20.259, 0),  # eE too large
    )

    # iterable with tuples (eNu, cosT, dSigma_dE(eNu, cosT))
    test_dSigma_dCosT_values = (
        (3, -0.50, 0.0),
        (3, 0.00, 0.0),
        (3, 0.50, 4.907277090752708e-22),
        (20, -0.99, 0.0),
        (20, 0.00, 0.0),
        (20, 0.99, 0.0), 
        (50, -0.50, 0.0),
        (50, 0.00, 0.0),
        (50, 0.50, 1.3389448150874471e-19),
    )

    # iterable with tuples (eNu, cosT, get_eE(eNu, cosT))
    test_get_eE_values = (
        (3, -0.50, 938.2768460503567),
        (3, 0.00, 938.27205),
        (3, 0.50, 938.2768460503567),
        (20, -0.99, 939.1077138141358),
        (20, 0.00, 938.27205),
        (20, 0.99, 939.1077138141358),
        (50, -0.50, 939.6042862101696),
        (50, 0.00, 938.27205),
        (50, 0.50, 939.6042862101696),
    )

    # iterable with tuples (eNu, bounds_eE(eNu)[0], bounds_eE(eNu)[1])
    test_bounds_eE_values = (
        (3, 938.27205, 938.2911123030725),
        (20, 938.27205, 939.0898184315932),
        (50, 938.27205, 943.087743536198),
    )

    # value of bounds_eNu[0]
    test_bounds_eNu_minvalue = 0.0


# ensure that unittest doesn't run tests in the base class, via https://stackoverflow.com/a/22836015
del CrossSectionTest

if __name__ == "__main__":
    unittest.main()
