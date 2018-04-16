from interaction_channels import ibd
from scipy import integrate

def dSigma_dE(eE, eNu):
    # exchange order of arguments and convert from MeV^2 to cm^2
    return ibd.dSigma_dE(eNu, eE) / (5.067731E10)**2

def sigma(eNu):
    return integrate.quad(dSigma_dE, *ibd.bounds_eE(eNu), args=(eNu))

# values from table 1 of Strumia/Vissani (2003), arXiv:astro-ph/0302055.
sigma_th = {2.01: 0.00351e-41,
            2.25: 0.00735e-41,
            2.51: 0.0127e-41,
            2.80: 0.0202e-41,
            3.12: 0.0304e-41,
            3.48: 0.0440e-41,
            3.89: 0.0619e-41,
            4.33: 0.0854e-41,
            4.84: 0.116e-41,
            5.40: 0.155e-41,
            6.02: 0.205e-41,
            6.72: 0.269e-41,
            7.49: 0.349e-41,
            8.36: 0.451e-41,
            8.83: 0.511e-41,
            9.85: 0.654e-41,
            11.0: 0.832e-41,
            12.3: 1.05e-41,
            13.7: 1.33e-41,
            15.3: 1.67e-41,
            17.0: 2.09e-41,
            19.0: 2.61e-41,
            21.2: 3.24e-41,
            23.6: 4.01e-41,
            26.4: 4.95e-41,
            29.4: 6.08e-41,
            32.8: 7.44e-41,
            36.6: 9.08e-41,
            40.9: 11.0e-41,
            43.2: 12.1e-41,
            48.2: 14.7e-41,
            53.7: 17.6e-41,
            59.9: 21.0e-41,
            66.9: 25.0e-41,
            74.6: 29.6e-41,
            83.2: 34.8e-41,
            92.9: 40.7e-41,
            104: 47.3e-41,
            160: 81.0e-41}

for eNu in sorted(sigma_th):
    sigma_sv = sigma_th[eNu]
    sigma_calc = sigma(eNu)[0]
    delta = (sigma_calc - sigma_sv) / sigma_sv
    if abs(delta) > 0.005:
        print "Inaccuracy at eNu =", eNu
        print "Delta = %.2f percent" % (delta * 100)
