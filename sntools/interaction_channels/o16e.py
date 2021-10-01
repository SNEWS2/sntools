"""Implementation of nu_e + 16O -> X + e-

Based on arXiv:1807.02367 (calculations) and arXiv:1809.08398 (fit).

That paper only gives the total cross section sigma(eNu), not the differential
c.s. dSigma/dE (eNu, eE), but since we assume eE = eNu - eG MeV, we can write
the differential c.s. as sigma(eNu) * delta(eNu - eG - eE).
However, numpy doesn't implement a delta distribution and numpy's (numerical)
integration doesn't play nice with sympy's (symbolic) DiracDelta, see:
https://stackoverflow.com/questions/36755487/diracdelta-not-giving-correct-result#36755974
Instead, below we implement an approximation to DiracDelta: a function that's
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""

from __future__ import division

from math import log10
import random

from sntools.event import Event

epsilon = 0.001  # for approximating DiracDelta distribution below

# Excitation energy and parameters a, b and c (Table 4 of arXiv:1809.08398)
fit_parameters = {1: [15.21, -40.008, 4.918, 1.036],
                  2: [22.47, -39.305, 4.343, 0.961],
                  3: [25.51, -39.655, 5.263, 1.236],
                  4: [29.35, -39.166, 3.947, 0.901]}


def generate_event(eNu, dirx, diry, dirz):
    """Return an event with the appropriate incoming/outgoing particles.

    Input:
        eNu: neutrino energy
        dirx, diry, dirz: direction of outgoing particle (normalized to 1)
    """
    eE = get_eE(eNu, dirz)

    evt = Event(1008016)
    evt.incoming_particles.append((12, eNu, 0, 0, 1))  # incoming neutrino
    evt.incoming_particles.append((8016, 14900, 0, 0, 1))  # oxygen nucleus at rest
    evt.outgoing_particles.append((11, eE, dirx, diry, dirz))  # outgoing electron
    return evt


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["e"]

# List with minimum & maximum energy of incoming neutrino.
bounds_eNu = [fit_parameters[1][0] + 0.8, 100]  # 0.8 MeV = Cherenkov threshold of electron


def bounds_eE(eNu, *args):
    """Return kinematic bounds for integration over eE.

    Input:
        eNu:  neutrino energy (in MeV)
        args: [ignore this]
    Output:
        list with minimum & maximum allowed energy of outgoing (detected) particle
    """
    # smallest eE is at largest (allowed) excitation energy
    for g in range(1, 5):
        if eNu > fit_parameters[g][0] + epsilon:
            eMin = eNu - fit_parameters[g][0] - epsilon

    # largest eE is at smallest excitation energy
    eMax = eNu - fit_parameters[1][0] + epsilon

    return [eMin, eMax]


def get_eE(eNu, cosT=0):
    """Return energy (in MeV) of outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (in MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    # find allowed excitation energies
    allowed = []
    for g in range(1, 5):
        if eNu > fit_parameters[g][0] + epsilon:
            eE = eNu - fit_parameters[g][0]
            sigma = partial_dSigma_dE(eNu, eE, g)
            allowed.append([eE, sigma])

    # choose from allowed eE with probability proportional to partial cross-section
    sigma_max = max([sigma for _, sigma in allowed])
    while True:
        eE, sigma = random.choice(allowed)
        if sigma > sigma_max * random.random():
            break
    return eE


def dSigma_dE(eNu, eE):
    """Return differential cross section in MeV^-2.

    Inputs:
        eNu: neutrino energy
        eE:  energy of outgoing (detected) particle
    """
    sigma = 0
    for g in range(1, 5):
        sigma += partial_dSigma_dE(eNu, eE, g)

    sigma *= (5.067731E10)**2  # convert cm^2 to MeV^-2, see http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    return sigma / (2 * epsilon)  # Ensure that integration over eE yields sigma


def partial_dSigma_dE(eNu, eE, g):  # eq. (4) of arXiv:1809.08398
    eG, a, b, c = fit_parameters[g]

    if abs(eNu - eE - eG) > epsilon:
        return 0

    d = log10(eNu**0.25 - eG**0.25)
    log_sigma = a + b * d + c * d**2
    return 10**log_sigma


def dSigma_dCosT(eNu, cosT):  # eq. (B7) of arXiv:hep-ph/0307050
    """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    if abs(cosT) > 1:
        return 0
    x = ((eNu - 15) / 25)**4
    return 1 - cosT * (1 + x) / (3 + x)


def _bounds_eNu(eE):
    """Min/max neutrino energy that can produce a given positron energy."""
    return (eE + fit_parameters[1][0] - epsilon, eE + fit_parameters[4][0] + epsilon)


def _opts(eNu, *args):
    """Options for numerical integration with `scipy.nquad`."""
    # values of eE where dSigma_dE(eNu, eE) has a discontinuity, to increase accuracy
    p = []
    for g in range(1, 5):
        if eNu > fit_parameters[g][0] + epsilon:
            p.append(eNu - fit_parameters[g][0] - epsilon)
            p.append(eNu - fit_parameters[g][0] + epsilon)

    return {'points': p}


def _opts2(eE, *args):
    """Options for numerical integration

    Values of eNu where dSigma_dE(eNu, eE) has a discontinuity, to increase accuracy
    """
    p = []
    for g in range(1, 5):
        p.append(eE + fit_parameters[g][0] - epsilon)
        p.append(eE + fit_parameters[g][0] + epsilon)

    return {'points': p}
