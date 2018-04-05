from math import pi, sqrt, log


targets_per_molecule = 2 # number of free protons per water molecule
pid = -11


'''
Particle physics section.
* cross section for neutrino-electron scattering
* directionality of scattered electron

Based on Strumia/Vissani (2003), arXiv:astro-ph/0302055.
'''
mN = 939.5654 # neutron mass (MeV)
mP = 938.2721 # proton mass (MeV)
mE = 0.5109989 # electron mass (MeV)
mPi = 139.57018 # pion mass (MeV)
delta = mN - mP
mAvg = (mP+mN) / 2
alpha = 1 / 137.036 # fine structure constant
gF = 1.16637e-11 # Fermi coupling constant
sigma0 = 2 * mP * gF**2 * 0.9746**2 / (8 * pi * mP**2) # from eqs. (3), (11)

def dSigma_dE(eNu, eE): # eqs. (11), (3)
    # above eq. (11)
    s_minus_u = 2*mP*(eNu+eE) - mE**2
    t = mN**2 - mP**2 - 2*mP*(eNu-eE)

    # eq. (7)
    x = t / (4*mAvg**2)
    y = 1 - t/710000
    z = 1 - t/1000000
    f1 = (1 - 4.706 * x) / ((1-x) * y**2)
    f2 = 3.706 / ((1-x) * y**2)
    g1 = -1.27 / z**2
    g2 = 2 * g1 * mAvg**2 / (mPi**2 - t)

    # Use NLO approximation (eq. (10)) for A, B and C. This is accurate to
    # better than 0.1% for eNu < 40 MeV (see line 4 in table 2).
    A = mAvg**2 * (f1**2 - g1**2) * (t - mE**2) \
        - mAvg**2 * delta**2 * (f1**2 + g1**2) \
        - 2 * mE**2 * mAvg * delta * g1 * (f1 + f2)
    B = t * g1 * (f1 + f2)
    C = (f1**2 + g1**2) / 4

    abs_M_squared = A - B * s_minus_u + C * s_minus_u**2 # eq. (5)
    rad_correction = alpha/pi * (6 + 3./2 * log(mP/(2*eE)) + 1.2 * (mE/eE)**1.5) # eq. (14)

    return sigma0 / eNu**2 * abs_M_squared * (1 + rad_correction)


# probability distribution for the angle at which the positron is emitted
def dSigma_dCosT(eNu, cosT): # eq. (20)
    epsilon = eNu / mP
    eE = get_eE(eNu, cosT)
    pE = sqrt(eE**2 - mE**2)
    dE_dCosT = pE * epsilon / (1 + epsilon * (1 - cosT * eE / pE))
    return dE_dCosT * dSigma_dE(eNu, eE)


def get_eE(eNu, cosT): # eq. (21)
    epsilon = eNu / mP
    kappa = (1 + epsilon)**2 - (epsilon * cosT)**2
    return ((eNu - delta_cm) * (1 + epsilon) + epsilon * cosT * sqrt((eNu - delta_cm)**2 - mE**2 * kappa)) / kappa


# Bounds for integration over eE
delta_cm = (mN**2 - mP**2 - mE**2) / (2*mP)
def bounds_eE(eNu, *args): # ignore additional arguments handed over by scipy.integrate.nquad()
    s = 2*mP*eNu + mP**2
    pE_cm = sqrt((s-(mN-mE)**2) * (s-(mN+mE)**2)) / (2*sqrt(s))
    eE_cm = (s-mN**2+mE**2) / (2*sqrt(s))

    eE_min = eNu - delta_cm - eNu/sqrt(s) * (eE_cm + pE_cm)
    eE_max = eNu - delta_cm - eNu/sqrt(s) * (eE_cm - pE_cm)
    return [eE_min, eE_max]

# Bounds for integration over eNu
eThr = ((mN+mE)**2 - mP**2) / (2*mP) # threshold energy for IBD: ca. 1.8 MeV
bounds_eNu = [eThr, 100]
