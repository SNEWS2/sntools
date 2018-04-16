from interaction_channels import ibd
from scipy import integrate

mev2cm = 1 / 5.067731E10

def sigma(eNu):
    func = lambda _eE: ibd.dSigma_dE(eNu, _eE) * mev2cm**2
    return integrate.quad(func, *ibd.bounds_eE(eNu))[0]

def eE_avg(eNu):
    func = lambda _eE: _eE * ibd.dSigma_dE(eNu, _eE) * mev2cm**2
    numerator = integrate.quad(func, *ibd.bounds_eE(eNu))[0]
    return numerator / sigma(eNu)

def cosT_avg(eNu):
    func = lambda _cosT: _cosT * ibd.dSigma_dCosT(eNu, _cosT) * mev2cm**2
    numerator = integrate.quad(func, -1, 1)[0]
    return numerator / sigma(eNu)


# values from table 1 of Strumia/Vissani (2003), arXiv:astro-ph/0302055.
sigma_th = {2.01: [0.00351e-41, 0.719, -0.021],
            2.25: [0.00735e-41, 0.952, -0.025],
            2.51: [0.0127e-41, 1.21, -0.027],
            2.80: [0.0202e-41, 1.50, -0.027],
            3.12: [0.0304e-41, 1.82, -0.027],
            3.48: [0.0440e-41, 2.18, -0.027],
            3.89: [0.0619e-41, 2.58, -0.026],
            4.33: [0.0854e-41, 3.03, -0.025],
            4.84: [0.116e-41, 3.52, -0.024],
            5.40: [0.155e-41, 4.08, -0.023],
            6.02: [0.205e-41, 4.69, -0.022],
            6.72: [0.269e-41, 5.38, -0.020],
            7.49: [0.349e-41, 6.15, -0.018],
            8.36: [0.451e-41, 7.00, -0.016],
            8.83: [0.511e-41, 7.46, -0.015],
            9.85: [0.654e-41, 8.47, -0.013],
            11.0: [0.832e-41, 9.58, -0.010],
            12.3: [1.05e-41, 10.8, -0.007],
            13.7: [1.33e-41, 12.2, -0.003],
            15.3: [1.67e-41, 13.7, 0.0006],
            17.0: [2.09e-41, 15.5, 0.005],
            19.0: [2.61e-41, 17.4, 0.010],
            21.2: [3.24e-41, 19.5, 0.015],
            23.6: [4.01e-41, 21.8, 0.021],
            26.4: [4.95e-41, 24.4, 0.028],
            29.4: [6.08e-41, 27.3, 0.036],
            32.8: [7.44e-41, 30.5, 0.044],
            36.6: [9.08e-41, 34.1, 0.054],
            40.9: [11.0e-41, 38.0, 0.065],
            43.2: [12.1e-41, 40.2, 0.070],
            48.2: [14.7e-41, 44.8, 0.083],
            53.7: [17.6e-41, 49.9, 0.097],
            59.9: [21.0e-41, 55.6, 0.113],
            66.9: [25.0e-41, 61.8, 0.131],
            74.6: [29.6e-41, 68.8, 0.151],
            83.2: [34.8e-41, 76.5, 0.173],
            92.9: [40.7e-41, 85.0, 0.198],
            104: [47.3e-41, 94.5, 0.225],
            160: [81.0e-41, 144, 0.361]}

for eNu in sorted(sigma_th):
    print "********** eNu =", eNu, "**********"
    sigma_sv = sigma_th[eNu][0]
    sigma_calc = sigma(eNu)
    delta = (sigma_calc - sigma_sv) / sigma_sv
    if abs(delta) > 0.005: # 0.5% = max. rounding error for 3 significant digits
        print "sigma: Delta = %.2f percent" % (delta * 100)

    eE_sv = sigma_th[eNu][1]
    eE_calc = eE_avg(eNu)
    delta = (eE_calc - eE_sv) / eE_sv
    if abs(delta) > 0.005:
        print "eE: Delta = %.2f percent" % (delta * 100)

    cosT_sv = sigma_th[eNu][2]
    cosT_calc = cosT_avg(eNu)
    delta = (cosT_calc - cosT_sv) / cosT_sv
    if abs(delta) > 0.01:
        print "cosT: Delta = %.2f percent" % (delta * 100)
