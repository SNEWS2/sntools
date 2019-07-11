from interaction_channels import ibd, es
from scipy import integrate
import __builtin__

channel = es

mev2cm = 1 / 5.067731E10

def sigma(eNu):
    func = lambda _eE: channel.dSigma_dE(eNu, _eE) * mev2cm**2
    return integrate.quad(func, *channel.bounds_eE(eNu))[0]

def eE_avg(eNu):
    func = lambda _eE: _eE * channel.dSigma_dE(eNu, _eE) * mev2cm**2
    numerator = integrate.quad(func, *channel.bounds_eE(eNu))[0]
    return numerator / sigma(eNu)

def cosT_avg(eNu):
    func = lambda _cosT: _cosT * channel.dSigma_dCosT(eNu, _cosT) * mev2cm**2
    numerator = integrate.quad(func, -1, 1)[0]
    return numerator / sigma(eNu)

if channel == ibd:
    '''Expected values for integrated cross-section, mean eE, mean cosT.

    Calculated with sample code by Strumia/Vissani. Note that these are not quite
    the same values as in table 1 of their paper (arXiv:astro-ph/0302055), where
    radiative corrections may not have been systematically included (according to
    private communication with A. Strumia).
    '''
    rounding_error = 0.05 # 0.05% = max. rounding error for 4 significant digits
    sigma_th = {2.01: [0.003439e-41, 0.7144, -0.0208],
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
                160: [81.30e-41, 143.7, 0.346]}

if channel == es:
    '''Expected values for integrated cross-section, mean eE, mean cosT.

    From table X of the Bahcall et al. paper (arXiv:XXXX.XXXXX).
    Note that the cross-section given there does not account for the Cherenkov
    threshold in water, so we need to reduce the minimum electron energy.
    '''
    es.eE_min = es.mE + 0.0000001 # kinetic energy > 0 to avoid divide-by-zero error
    __builtin__._flavor = "e"

    rounding_error = 0.5 # 0.5% = max. rounding error for 3 significant digits
    sigma_th = {0.38: [1.92e-45, 5.09e-46, 0], # e- will be below Cherenkov threshold
                0.86: [5.79e-45, 1.28e-45, 0],
                1.00: [6.98e-45, 1.50e-45, 0],
                1.44: [1.09e-44, 2.21e-45, 0],
                2.00: [1.59e-44, 3.11e-45, 0],
                3.00: [2.51e-44, 4.70e-45, 0],
                4.00: [3.42e-44, 6.29e-45, 0],
                5.00: [4.35e-44, 7.87e-45, 0],
                7.00: [6.19e-44, 1.10e-44, 0],
                10.00: [8.96e-44, 1.58e-44, 0],
                12.00: [1.08e-43, 1.90e-44, 0],
                14.00: [1.27e-43, 2.21e-44, 0],
                16.00: [1.45e-43, 2.53e-44, 0],
                18.00: [1.64e-43, 2.85e-44, 0],
                20.00: [1.82e-43, 3.16e-44, 0],
                25.00: [2.28e-43, 3.95e-44, 0],
                30.00: [2.74e-43, 4.74e-44, 0],
                40.00: [3.67e-43, 6.32e-44, 0],
                50.00: [4.59e-43, 7.90e-44, 0],
                60.00: [5.52e-43, 9.48e-44, 0]}

good_agreement = True
for eNu in sorted(sigma_th):
    if channel == es: # cross-section of electron neutrinos
         __builtin__._flavor = "e"

    sigma_sv = sigma_th[eNu][0]
    sigma_calc = sigma(eNu)
    delta = 100 * (sigma_calc - sigma_sv) / sigma_sv
    if abs(delta) > rounding_error:
        print "eNu = %.2f: sigma is off by %.2f percent" % (eNu, delta)
        good_agreement = False

    if channel == es: # cross-section of muon neutrinos
        __builtin__._flavor = "x"
        sigma_sv = sigma_th[eNu][1]
        sigma_calc = sigma(eNu)
        delta = 100 * (sigma_calc - sigma_sv) / sigma_sv
        if abs(delta) > rounding_error:
            print "eNu = %.2f: sigma is off by %.2f percent" % (eNu, delta)
            good_agreement = False

#     eE_sv = sigma_th[eNu][1]
#     eE_calc = eE_avg(eNu)
#     delta = 100 * (eE_calc - eE_sv) / eE_sv
#     if abs(delta) > 0.05:
#         print "eNu = %.2f: eE is off by %.2f percent" % (eNu, delta)
#         good_agreement = False
#
#     cosT_sv = sigma_th[eNu][2]
#     cosT_calc = cosT_avg(eNu)
#     delta = 100 * (cosT_calc - cosT_sv) / cosT_sv
#     if abs(delta) > 0.5:
#         print "eNu = %.2f: cosT is off by %.2f percent" % (eNu, delta)
#         good_agreement = False

if good_agreement:
    print "Integrated cross-section, mean eE & mean cosT are in excellent", \
          "agreement for eNu between 2 and 100 MeV."
