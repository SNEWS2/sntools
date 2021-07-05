from __future__ import print_function
import unittest
from scipy import integrate

from sntools.interaction_channels import ibd


mev2cm = 1 / 5.067731e10

c = ibd.Channel('eb')


def sigma(eNu):
    func = lambda _eE: c.dSigma_dE(eNu, _eE) * mev2cm ** 2
    return integrate.quad(func, *c.bounds_eE(eNu))[0]


def eE_avg(eNu):
    func = lambda _eE: _eE * c.dSigma_dE(eNu, _eE) * mev2cm ** 2
    numerator = integrate.quad(func, *c.bounds_eE(eNu))[0]
    return numerator / sigma(eNu)


def cosT_avg(eNu):
    func = lambda _cosT: _cosT * c.dSigma_dCosT(eNu, _cosT) * mev2cm ** 2
    numerator = integrate.quad(func, -1, 1)[0]
    return numerator / sigma(eNu)


"""Expected values for integrated cross-section, mean eE, mean cosT.

Calculated with sample code by Strumia/Vissani. Note that these are not quite
the same values as in table 1 of their paper (arXiv:astro-ph/0302055), where
radiative corrections may not have been systematically included (according to
private communication with A. Strumia).
"""
sigma_th = {
    2.01: [0.003439e-41, 0.7144, -0.0208],
    2.25: [0.007381e-41, 0.9536, -0.0254],
    2.51: [0.01279e-41, 1.213, -0.0270],
    2.80: [0.02027e-41, 1.501, -0.0274],
    3.12: [0.03032e-41, 1.820, -0.0273],
    3.48: [0.04385e-41, 2.178, -0.0269],
    3.89: [0.06213e-41, 2.585, -0.0262],
    4.33: [0.08510e-41, 3.022, -0.0253],
    4.84: [0.1160e-41, 3.527, -0.0242],
    5.40: [0.1552e-41, 4.082, -0.0230],
    6.02: [0.2050e-41, 4.695, -0.0216],
    6.72: [0.2690e-41, 5.387, -0.0200],
    7.49: [0.3489e-41, 6.146, -0.0181],
    8.36: [0.4509e-41, 7.003, -0.0161],
    8.83: [0.5111e-41, 7.465, -0.0150],
    9.85: [0.6538e-41, 8.466, -0.0125],
    11.0: [0.8340e-41, 9.593, -0.00975],
    12.3: [1.062e-41, 10.86, -0.00662],
    13.7: [1.334e-41, 12.23, -0.00323],
    15.3: [1.680e-41, 13.78, 0.000664],
    17.0: [2.085e-41, 15.43, 0.00481],
    19.0: [2.608e-41, 17.36, 0.00971],
    21.2: [3.241e-41, 19.48, 0.0151],
    23.6: [3.996e-41, 21.77, 0.0210],
    26.4: [4.954e-41, 24.44, 0.0280],
    29.4: [6.068e-41, 27.29, 0.0354],
    32.8: [7.430e-41, 30.50, 0.0439],
    36.6: [9.066e-41, 34.06, 0.0534],
    40.9: [11.05e-41, 38.08, 0.0642],
    43.2: [12.16e-41, 40.21, 0.0700],
    48.2: [14.67e-41, 44.83, 0.0826],
    53.7: [17.59e-41, 49.88, 0.0965],
    59.9: [21.03e-41, 55.53, 0.112],
    66.9: [25.06e-41, 61.86, 0.130],
    74.6: [29.63e-41, 68.78, 0.149],
    83.2: [34.85e-41, 76.46, 0.171],
    92.9: [40.81e-41, 85.06, 0.194],
    104: [47.68e-41, 94.84, 0.221],
    160: [81.30e-41, 143.7, 0.346],
}


class ExtendedIBDTest(unittest.TestCase):
    def test_sigma(self):
        for eNu in sorted(sigma_th):
            sigma_sv = sigma_th[eNu][0]
            # 0.0005 = max. rounding error for 4 significant digits
            self.assertAlmostEqual(sigma(eNu), sigma_sv, delta=0.0005 * sigma_sv)

    def test_eE(self):
        for eNu in sorted(sigma_th):
            eE_sv = sigma_th[eNu][1]
            self.assertAlmostEqual(eE_avg(eNu), eE_sv, delta=0.0005 * eE_sv)

    def test_cosT(self):
        for eNu in sorted(sigma_th):
            cosT_sv = sigma_th[eNu][2]
            self.assertAlmostEqual(cosT_avg(eNu), cosT_sv, delta=abs(0.005 * cosT_sv))
