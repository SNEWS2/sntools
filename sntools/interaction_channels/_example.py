'''
Sample implementation of an interaction channel.

Most of the actual work is done in `channel.py`. Here, you need to provide
constants and functions that characterize this interaction channel.
See the docstrings below for detailed descriptions.

If you need to define helper functions or constants, you can do so at the bottom
of this file, where some commonly used constants are already provided.
'''

from sntools.event import Event

'''
generate_event(eNu, direction):
Generate an event with appropriate incoming/outgoing particles.
Input:
    eNu: neutrino energy
    dirx, diry, dirz: direction of outgoing particle (normalized to 1)
Output:
    event object
'''
def generate_event(eNu, dirx, diry, dirz):
    eE = get_eE(eNu, dirz)

    # replace `0` with interaction channel code ±CZZZAAA, where
    # + is for neutrinos, - for antineutrinos
    # C is 1 for charged-current and 2 for neutral-current interaction
    # AAAZZZ depends on the nucleus, e.g. 006012 for carbon-12 or 008016 for oxygen-16
    evt = Event(0)
#     evt.incoming_particles.append([12, eNu, 0, 0, 1]) # incoming neutrino
#     evt.incoming_particles.append([11, mE, 0, 0, 1]) # electron at rest
#     evt.outgoing_particles.append([11, eE, dirx, diry, dirz]) # outgoing electron
    return evt


'''
possible_flavors:
which neutrino flavors ("e", "eb", "x", "xb") interact in this channel
'''
possible_flavors = ["e"]


'''
dSigma_dE(eNu, eE):
Differential cross section.
Input:
    eNu: neutrino energy
    eE:  energy of outgoing (detected) particle
Output:
    one floating point number
'''
def dSigma_dE(eNu, eE):
    return None


'''
dSigma_dCosT(eNu, cosT):
Distribution of the angle at which the outgoing (detected) particle is emitted.
Input:
    eNu:  neutrino energy (MeV)
    cosT: cosine of the angle between neutrino and outgoing (detected) particle
Output:
    one floating point number
'''
def dSigma_dCosT(eNu, cosT):
    return None
#     eE = get_eE(eNu, cosT)
#     dE_dCosT = None
#     return dE_dCosT * dSigma_dE(eNu, eE)


'''
get_eE(eNu, cosT):
Energy of outgoing (detected) particle.
Input:
    eNu:  neutrino energy (MeV)
    cosT: cosine of the angle between neutrino and outgoing (detected) particle
Output:
    one floating point number
'''
def get_eE(eNu, cosT):
    return None


'''
bounds_eE(eNu, *args):
Kinematical bounds for integration over eE.
Input:
    eNu:  neutrino energy (MeV)
    args: [implementation detail; ignore this]
Output:
    list with minimum & maximum allowed energy of outgoing (detected) particle
'''
def bounds_eE(eNu, *args):
    return [None, None]
#     return [eE_min(eNu), eE_max(eNu)]
#
# def eE_min(eNu):
#     return None
#
# def eE_max(eNu):
#     return None


'''
bounds_eNu
List with minimum & maximum energy of incoming neutrino. The minimum energy is
given by the threshold energy for the interaction (where applicable),
while the maximum energy is given by the supernova neutrino flux.
'''
bounds_eNu = [None, 100]
# e_threshold = 0 # threshold energy for current channel
# bounds_eNu = [e_threshold, 100]

# minimum/maximum neutrino energy that can produce a given positron energy
# Optional. Can reduce numerical inaccuracy when integrating over eNu.
def _bounds_eNu(eE):
    return bounds_eNu


'''
End of required values.

If you need to define helper functions or constants, you can do so below.
Some commonly needed constants are already provided.
'''
# mN = 939.5654 # neutron mass (MeV)
# mP = 938.2721 # proton mass (MeV)
# mE = 0.5109989 # electron mass (MeV)
# alpha = 1 / 137.036 # fine structure constant
# gF = 1.16637e-11 # Fermi coupling constant
