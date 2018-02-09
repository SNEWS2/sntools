'''
Implementation of nu_e + 16O -> X + e-

Based on Appendix B.3 of hep-ph/0307050.
That paper only gives the toal cross section sigma(eNu), not the differential
c.s. dSigma/dE (eNu, eE), but since we assume eE = eNu - 15 MeV, we can write
the differential c.s. as sigma(eNu) * delta(eNu - 15 - eE).
However, numpy doesn't implement a delta distribution and numpy's (numerical)
integration doesn't play nice with sympy's (symbolic) DiracDelta, see:
https://stackoverflow.com/questions/36755487/diracdelta-not-giving-correct-result#36755974
Instead, below we implement an approximation to DiracDelta: a function that's
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
'''

e_thr = 15 # energy threshold for this reaction (~nuclear binding energy in 16O)
epsilon = 0.001 # for approximating DiracDelta distribution below


'''
targets_per_molecule:
number of interaction targets per water molecule
(i.e. 2 free protons, 1 oxygen nucleus or 10 electrons)
'''
targets_per_molecule = 1


'''
pid:
ID of the outgoing (detected) particle, using Particle Data Group conventions
(e.g. electron = 11, positron = -11)
'''
pid = 11


'''
bounds_eNu
List with minimum & maximum energy of incoming neutrino. The minimum energy is
typically given by the threshold energy for the interaction, while the maximum
energy is given by the supernova neutrino flux.
'''
bounds_eNu = [e_thr + 0.8, 100] # 0.8 MeV = Cherenkov threshold of electron


'''
bounds_eE(eNu, *args):
Kinematical bounds for integration over eE.
Input:
    eNu:  neutrino energy (MeV)
    args: [ignore this]
Output:
    list with minimum & maximum allowed energy of outgoing (detected) particle
'''
def bounds_eE(eNu, *args):
    return [get_eE(eNu) - epsilon, get_eE(eNu) + epsilon]


'''
get_eE(eNu, cosT):
Energy of outgoing (detected particle).
Input:
    eNu:  neutrino energy (MeV)
    cosT: cosine of the angle between neutrino and outgoing (detected) particle
Output:
    one floating point number
'''
def get_eE(eNu, cosT=0):
    return eNu - e_thr


'''
dSigma_dE(eNu, eE):
Differential cross section.
Input:
    eNu: neutrino energy
    eE:  energy of outgoing (detected) particle
Output:
    one floating point number
'''
def dSigma_dE(eNu, eE): # eq. (B6)
    if abs(get_eE(eNu) - eE) > epsilon:
        # This should never be called since we set bounds_eE() accordingly above
        # ... but just in case:
        return 0

    sigma0 = 4.7E-40 * (5.067731E10)**2 # convert cm^2 to MeV^-2, see http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    sigma = sigma0 * (eNu**0.25 - 15**0.25)**6
    return sigma / (2*epsilon) # Ensure that integration over eE yields sigma


'''
dSigma_dCosT(eNu, cosT):
Distribution of the angle at which the outgoing (detected) particle is emitted.
Input:
    eNu:  neutrino energy (MeV)
    cosT: cosine of the angle between neutrino and outgoing (detected) particle
Output:
    one floating point number
'''
def dSigma_dCosT(eNu, cosT): # eq. (B7)
    x = (get_eE(eNu, cosT) / 25)**4
    return 1 - cosT * (1+x)/(3+x)
