from math import pi, sqrt, log
from scipy import integrate

# Note: `es.py` uses `__builtin__._flavor` (which is set in `genevts.py`)

targets_per_molecule = 10 # number of electrons per water molecule
pid = 11


'''
Particle physics section.
* cross section for neutrino-electron scattering
* directionality of scattered electron

Based on the Appendices of Bahcall et al. 1995 (https://doi.org/10.1103/PhysRevD.51.6146).
This calculation includes radiative corrections from QCD & QED effects.

For differences between neutrinos/antineutrinos and a derivation of directionality,
see https://www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
Note that it uses different conventions (e.g. minus signs) from Bahcall et al.!
'''
sin2theta_w = 0.2317 # weak mixing angle
alpha = 1 / 137.036 # fine structure constant
mE = 0.5109989 # electron mass (MeV)
gF = 1.16637e-11 # Fermi coupling constant
rho_NC = 1.0126 # numerical factor from Bahcall et al.

def dSigma_dE(eNu, eE):
    return 2*mE*gF**2 / pi * (
               gL(eE)**2 * (1 + alpha/pi * f1(eNu, eE))
               + gR(eE)**2 * ((1-z(eNu, eE))**2 + alpha/pi * f2(eNu, eE))
               - gR(eE) * gL(eE) * mE/eNu * z(eNu, eE) * (1 + alpha/pi * f3(eNu, eE))
               )

# Appendix A: Radiative Corrections
def l(eE):
    return sqrt(eE**2 - mE**2)
def beta(eNu, eE):
    return l(eE)/eNu
def T(eE): # kinetic energy of recoil electron
    return eE - mE
def z(eNu, eE):
    return T(eE)/eNu
def x(eE):
    return sqrt(1 + 2*mE/T(eE))
def I(eE):
    return 1./6 * (1./3 + (3 - x(eE)**2) * (x(eE)/2. * log((x(eE)+1)/(x(eE)-1)) - 1))
def k(eE):
    if _flavor in ("e", "eb"):
        return 0.9791 + 0.0097 * I(eE) #+/-0.0025
    if _flavor in ("x", "xb"):
        return 0.9970 - 0.00037 * I(eE) #+/-0.0025

def g1(eE):
    return rho_NC * (0.5 - k(eE) * sin2theta_w)
def g2(eE):
    return -rho_NC * k(eE) * sin2theta_w

def gL(eE):
    if   _flavor == "e":  return g1(eE) - 1
    elif _flavor == "eb": return g2(eE)
    elif _flavor == "x":  return g1(eE)
    elif _flavor == "xb": return g2(eE)

def gR(eE):
    if   _flavor == "e":  return g2(eE)
    elif _flavor == "eb": return g1(eE) - 1
    elif _flavor == "x":  return g2(eE)
    elif _flavor == "xb": return g1(eE)

# Appendix B: QED Effects
def spence(n):
    return integrate.quad(lambda t: log(abs(1-t))/t, 0, n) [0]

def f1(eNu, eE): # fMinus(z)
    _z = z(eNu, eE)
    return (eE/l(eE) * log((eE+l(eE))/mE) - 1) \
                * (2*log(1 - _z - mE/(eE+l(eE))) - log(1 - _z) - log(_z)/2. - 5./12) \
           + 0.5 * (spence(_z) - spence(beta(eNu, eE))) \
           - 0.5 * log(1-_z)**2 - (11./12 + _z/2.) * log(1-_z) \
           + _z * (log(_z) + 0.5 * log(2*eNu / mE)) \
           - (31./18 + 1./12 * log(_z)) * beta(eNu, eE) \
           - 11./12 * _z + _z**2 / 24.

def f2(eNu, eE): # (1-z)**2 * fPlus(z)
    _z = z(eNu, eE)
    return (eE/l(eE) * log((eE+l(eE))/mE) - 1) \
                * ((1-_z)**2 * (2*log(1-_z-mE/(eE+l(eE))) - log(1-_z) - log(_z)/2. - 2./3) - (_z**2 * log(_z) + 1 - _z)/2.) \
           - (1-_z)**2 / 2. * (log(1-_z)**2 + beta(eNu, eE) * (spence(1-_z) - log(_z)*log(1-_z))) \
           + log(1-_z) * (_z**2 / 2. * log(_z) + (1-_z)/3. * (2*_z - 0.5)) \
           - _z**2 / 2. * spence(1-_z) - _z * (1-2*_z)/3 * log(_z) - _z * (1- _z)/6 \
           - beta(eNu, eE)/12. * (log(_z) + (1-_z) * (115 - 109 * _z)/6.)

def f3(eNu, eE): # fPlusMinus(z)
    return (eE/l(eE) * log((eE+l(eE))/mE) - 1) * 2 * log(1 - z(eNu, eE) - mE/(eE + l(eE)))

# energy of electron scattered into direction cosT by a neutrino with energy eNu
def get_eE(eNu, cosT):
    return mE + (2 * mE * eNu**2 * cosT**2) / ((mE + eNu)**2 - eNu**2 * cosT**2)

# distribution of scattering angles
def dSigma_dCosT(eNu, cosT):
    if cosT < 0: # backward scattering is kinematically impossible
        return 0

    dE_dCosT = 4 * mE * eNu**2 * (mE+eNu)**2 * cosT / ((mE+eNu)**2 - eNu**2 * cosT**2)**2
    eE = get_eE(eNu, cosT)
    return dE_dCosT * dSigma_dE(eNu, eE)


# Bounds for integration over eE
eE_min = 0.77 # Cherenkov threshold in water (refraction index n=1.34)
def eE_max(eNu):
    return mE + 2*eNu**2 / (2*eNu + mE) # this is get_eE(eNu, cosT=1)
def bounds_eE(eNu, *args): # ignore additional arguments handed over by integrate.nquad()
    return [eE_min, eE_max(eNu)]


# Bounds for integration over eNu
def eNu_min(eE):
    return T(eE)/2. * (1 + sqrt(1 + 2*mE/T(eE))) # inversion of eE_max(eNu)
eNu_max = 50
bounds_eNu = [eNu_min(eE_min), eNu_max]
