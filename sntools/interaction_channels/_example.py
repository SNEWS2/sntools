"""Sample implementation of an interaction channel.

Most of the actual work is done in `channel.py`. Here, you need to provide
constants and functions that characterize this interaction channel.
See the comments below for detailed descriptions.
"""

from sntools.event import Event


def generate_event(eNu, dirx, diry, dirz):
    """Return an event with the appropriate incoming/outgoing particles.

    Input:
        eNu: neutrino energy
        dirx, diry, dirz: direction of outgoing particle (normalized to 1)
    """
    eE = get_eE(eNu, dirz)

    # replace `0` with interaction channel code ±CZZZAAA, where
    # + is for neutrinos, - for antineutrinos
    # C is 1 for charged-current and 2 for neutral-current interaction
    # AAAZZZ depends on the nucleus, e.g. 006012 for carbon-12 or 008016 for oxygen-16
    evt = Event(0)
#     evt.incoming_particles.append([12, eNu, 0, 0, 1])  # incoming neutrino
#     evt.incoming_particles.append([11, mE, 0, 0, 1])  # electron at rest
    evt.outgoing_particles.append([11, eE, dirx, diry, dirz])  # outgoing electron
    return evt


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["e"]


def dSigma_dE(eNu, eE):
    """Return differential cross section in MeV^-2.

    Inputs:
        eNu: neutrino energy
        eE:  energy of outgoing (detected) particle
    """
    sigma = 0
    sigma *= (5.067731E10)**2  # if necessary, convert cm^2 to MeV^-2, see http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    return sigma


def dSigma_dCosT(eNu, cosT):
    """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    return 0.0


def get_eE(eNu, cosT):
    """Return energy (in MeV) of outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (in MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    return eNu


def bounds_eE(eNu, *args):
    """Return kinematic bounds for integration over eE.

    Input:
        eNu:  neutrino energy (in MeV)
        args: [ignore this]
    Output:
        list with minimum & maximum allowed energy of outgoing (detected) particle
    """
    return [0, eNu]


# List with minimum & maximum energy of incoming neutrino.
# The minimum is usually given by the energy threshold for the interaction,
# while the maximum energy is given by the supernova neutrino flux.
bounds_eNu = [0, 100]


# End of required values.
# The following functions are optional but may improve numerical accuracy.

def _bounds_eNu(eE):
    """Min/max neutrino energy that can produce a given positron energy.

    Optional. Can reduce numerical inaccuracy when integrating over eNu.
    """
    return bounds_eNu


def _opts(eNu, *args):
    """Options for numerical integration with `scipy.nquad`.

    Optional. Return values of eE where dSigma_dE(eNu, eE) has a discontinuity, to improve numerical inaccuracy.
    """
    return {'points': []}
