"""Implementation of anti-nu_e + 16O -> X + e+

Very schematic implementation. Based on data in PRD 66,013007 and
somewhat analogous to Appendix B.3 of hep-ph/0307050.
That paper only gives the total cross section sigma(eNu), not the differential
c.s. dSigma/dE (eNu, eE), but since we assume eE = eNu - 11.4 MeV, we can write
the differential c.s. as sigma(eNu) * delta(eNu - 11.4 - eE).
However, numpy doesn't implement a delta distribution and numpy's (numerical)
integration doesn't play nice with sympy's (symbolic) DiracDelta, see:
https://stackoverflow.com/questions/36755487/diracdelta-not-giving-correct-result#36755974
Instead, below we implement an approximation to DiracDelta: a function that's
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""
from sntools.event import Event


def generate_event(eNu, dirx, diry, dirz):
    """Return an event with the appropriate incoming/outgoing particles.

    Input:
        eNu: neutrino energy
        dirx, diry, dirz: direction of outgoing particle (normalized to 1)
    """
    eE = get_eE(eNu, dirz)

    evt = Event(-1008016)
    evt.incoming_particles.append((-12, eNu, 0, 0, 1))  # incoming neutrino
    evt.incoming_particles.append((8016, 14900, 0, 0, 1))  # oxygen nucleus at rest
    evt.outgoing_particles.append((-11, eE, dirx, diry, dirz))  # outgoing positron
    return evt


e_thr = 11.4  # energy threshold for this reaction
epsilon = 0.001  # for approximating DiracDelta distribution below


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["eb"]


# List with minimum & maximum energy of incoming neutrino.
bounds_eNu = [e_thr + 0.8, 100]  # 0.8 MeV = Cherenkov threshold of electron


def bounds_eE(eNu, *args):
    """Return kinematic bounds for integration over eE.

    Input:
        eNu:  neutrino energy (in MeV)
        args: [ignore this]
    Output:
        list with minimum & maximum allowed energy of outgoing (detected) particle
    """
    return [get_eE(eNu) - epsilon, get_eE(eNu) + epsilon]


def get_eE(eNu, cosT=0):
    """Return energy (in MeV) of outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (in MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    return eNu - e_thr


def dSigma_dE(eNu, eE):
    """Return differential cross section in MeV^-2.

    Inputs:
        eNu: neutrino energy
        eE:  energy of outgoing (detected) particle
    """
    if abs(get_eE(eNu) - eE) > epsilon:
        # This should never be called since we set bounds_eE() accordingly above
        # ... but just in case:
        return 0

    # parameters from my own fit to data in PRD 66,013007 (Table 1)
    if eNu < 54:
        sigma0 = 9.35E-44
        a = 0.497
        b = 4.753
    else:
        sigma0 = 7.91E-41
        a = 0.2616
        b = 5.18
    # see http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    cm2mev = 5.067731E10
    sigma = sigma0 * (eNu**a - e_thr**a)**b * cm2mev**2
    return sigma / (2 * epsilon)  # Ensure that integration over eE yields sigma


def dSigma_dCosT(eNu, cosT):
    """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    # Plots in PRD 36,2283 show this behaves roughly similar to the analogous
    # nu_e reaction, so we use the same approximation. (hep-ph/0307050, eq. B7)
    x = (get_eE(eNu, cosT) / 25)**4
    return 1 - cosT * (1 + x) / (3 + x)


def _bounds_eNu(eE):
    """Min/max neutrino energy that can produce a given positron energy."""
    return (eE + e_thr - epsilon, eE + e_thr + epsilon)
