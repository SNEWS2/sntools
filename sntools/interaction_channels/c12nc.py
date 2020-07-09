'''
Implementation of nu + 12C -> nu' + 12C*, 12C* -> 12C + \gamma

Based on Donnelly & Peccei, Phys.Rept. 50 (1979) 1, eq. (4.53). We assume
\beta * \kappa = 1.11, based on an experimental determination, see
Armbruster et al. (KARMEN Collaboration), Phys.Lett.B423 (1998), p. 15.
This is higher than the SM prediction (=1) but consistent with SNOwGLoBES.

To determine the the differential cross section dSigma/dE (eNu, eE) from the
total cross section, we approximate a DiracDelta function with one that is
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
'''

from sntools.event import Event

e_thr = 15.11 # energy threshold of this reaction
epsilon = 0.001 # for approximating DiracDelta distribution below


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
    # Note: `builtins._flavor` is set in `genevts.py`
    nu_flv = {'e':12, 'eb':-12, 'x':14, 'xb':-14}[_flavor]

    evt = Event(2006012 if nu_flv > 0 else -2006012)
    evt.incoming_particles.append([nu_flv, eNu, 0, 0, 1]) # incoming nu
    evt.incoming_particles.append((6012, 11178, 0, 0, 1)) # carbon nucleus at rest
    evt.outgoing_particles.append([22, e_thr, dirx, diry, dirz]) # emitted gamma
#     evt.outgoing_particles.append([nu_flv, eNu-e_thr, 0, 0, 1]) # outgoing nu
    return evt


'''
possible_flavors:
which neutrino flavors ("e", "eb", "x", "xb") interact in this channel
'''
possible_flavors = ("e", "eb", "x", "xb")


'''
bounds_eNu
List with minimum & maximum energy of incoming neutrino. The minimum energy is
typically given by the threshold energy for the interaction, while the maximum
energy is given by the supernova neutrino flux.
'''
bounds_eNu = (e_thr, 100)


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
Energy of outgoing (detected) particle.
Input:
    eNu:  neutrino energy (MeV)
    cosT: cosine of the angle between neutrino and outgoing (detected) particle
Output:
    one floating point number
'''
def get_eE(eNu, cosT=0):
    return e_thr


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
    if eNu < e_thr or abs(get_eE(eNu) - eE) > epsilon:
        # This should never happen, since we set bounds for eE and eNu accordingly above
        # ... but just in case:
        return 0

    sigma0 = 1.08E-38 * (5.067731E10)**2 # convert cm^2 to MeV^-2, see http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
    sigma = sigma0 * ((eNu-e_thr) / 939)**2
    sigma *= 1.11**2 # coupling constant \beta (see discussion above)
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
def dSigma_dCosT(eNu, cosT):
    # Energy dependence is unclear, so we use a constant value for now.
    if abs(cosT) > 1:
        return 0
    return 0.5


# min/max neutrino energy to produce outgoing particle with observed energy
def _bounds_eNu(eE):
    return bounds_eNu
