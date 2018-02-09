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
delta = mN-mP
mAvg = (mP+mN)/2
alpha = 1 / 137.036 # fine structure constant
gF = 1.16637e-11 # Fermi coupling constant
sigma0 = 2 * mP * gF**2 * 0.9746**2 / (8 * pi * mP**2) # from eqs. (3), (11)

def t_eNu_eE(eNu, eE): # above eq. (11)
    return mN**2 - mP**2 - 2*mP*(eNu-eE)
def x(eNu, eE):
    return t_eNu_eE(eNu, eE) / (4*mAvg**2)
def y(eNu, eE):
    return 1 - t_eNu_eE(eNu, eE)/710000
def z(eNu, eE):
    return 1 - t_eNu_eE(eNu, eE)/1000000
def f1(eNu, eE): # eq. (7)
    return (1 - 4.706*x(eNu, eE)) / ((1-x(eNu, eE)) * y(eNu, eE)**2)
def f2(eNu, eE): # eq. (7)
    return 3.706 / ((1-x(eNu, eE)) * y(eNu, eE)**2)
def g1(eNu, eE): # eq. (7)
    return -1.27 / z(eNu, eE)**2
def g2(eNu, eE): # eq. (7)
    return 2 * g1(eNu, eE) * mAvg**2 / (mPi**2 - t_eNu_eE(eNu, eE))

# Use NLO approximation (eq. (10)) for A, B and C. This is accurate to
# better than 0.1% for eNu < 40 MeV (see line 4 in table 2).
def AM(eNu, eE):
    return (mAvg**2 * (f1(eNu, eE)**2 - g1(eNu, eE)**2) * (t_eNu_eE(eNu, eE)-mE**2)) - (mAvg**2 * delta**2 * (f1(eNu, eE)**2 + g1(eNu, eE)**2)) - (2 * mE**2 * mAvg * delta * g1(eNu, eE) *(f1(eNu, eE)+f2(eNu, eE)))

def BM(eNu, eE):
    return t_eNu_eE(eNu, eE)*g1(eNu, eE)*(f1(eNu, eE)+f2(eNu, eE))

def CM(eNu, eE):
    return ((f1(eNu, eE)**2) + (g1(eNu, eE)**2))/4

def sMinusU(eNu, eE): # above eq. (11)
    return 2*mP*(eNu+eE) - mE**2

def absMsquared(eNu, eE): # eq. (5)
    return AM(eNu, eE) - sMinusU(eNu, eE) * BM(eNu, eE) + sMinusU(eNu, eE)**2 * CM(eNu, eE)

def rad_correction(eE): # eq. (14)
    return 1 + alpha/pi * (6 + 3./2 * log(mP / (2 * eE)) + 1.2 * (mE/eE)**1.5)

def dSigma_dE(eNu, eE): # eqs. (11), (3)
    return sigma0 / eNu**2 * absMsquared(eNu, eE) * rad_correction(eE)


# probability distribution for the angle at which the positron is emitted
def dSigma_dCosT(eNu, cosT): # eq. (20)
    epsilon = eNu/mP
    eE = get_eE(eNu, cosT)
    pE = sqrt(eE**2 - mE**2)
    dE_dCosT = pE * epsilon / (1 + epsilon * (1 - cosT * eE / pE))
    return dE_dCosT * dSigma_dE(eNu, eE)


def get_eE(eNu, cosT): # eq. (21)
    epsilon = eNu/mP
    kappa = (1 + epsilon)**2 - (epsilon * cosT)**2
    return ((eNu - delta_cm) * (1 + epsilon) + epsilon * cosT * sqrt((eNu - delta_cm)**2 - mE**2 * kappa)) / kappa


# calculate range for eE from eNu in center-of-mass (cm) frame
def s(eNu):
    return 2*mP*eNu + mP**2
def pE_cm(eNu):
    return (sqrt((s(eNu)-(mN-mE)**2)*(s(eNu)-(mN+mE)**2)))/(2*sqrt(s(eNu)))
def eE_cm(eNu):
    return (s(eNu)-mN**2+mE**2)/(2*sqrt(s(eNu)))

# Bounds for integration over eE
delta_cm = (mN**2 - mP**2 - mE**2) / (2*mP)
def eE_min(eNu): # eq. (12)
    return eNu - delta_cm - (eNu/sqrt(s(eNu)) * (eE_cm(eNu) + pE_cm(eNu)))
def eE_max(eNu): # eq. (12)
    return eNu - delta_cm - (eNu/sqrt(s(eNu)) *(eE_cm(eNu) - pE_cm(eNu)))
def bounds_eE(eNu, *args): # ignore additional arguments handed over by scipy.integrate.nquad()
    return [eE_min(eNu)+1, eE_max(eNu)+1]

# Bounds for integration over eNu
eThr = ((mN+mE)**2 - mP**2) / (2*mP) # threshold energy for IBD
bounds_eNu = [eThr, 100]
