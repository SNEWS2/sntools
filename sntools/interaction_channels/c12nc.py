"""Implementation of nu + 12C -> nu' + 12C*, 12C* -> 12C + gamma

Based on Donnelly & Peccei, Phys.Rept. 50 (1979) 1, eq. (4.53).
We assume beta * kappa = 1.11, based on an experimental determination, see
Armbruster et al. (KARMEN Collaboration), Phys.Lett.B423 (1998), p. 15.
This is higher than the SM prediction (=1) but consistent with SNOwGLoBES.

To determine the the differential cross section dSigma/dE (eNu, eE) from the
total cross section, we approximate a DiracDelta function with one that is
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""

from sntools.event import Event
from sntools.interaction_channels import BaseChannel

e_thr = 15.11  # energy threshold of this reaction
epsilon = 0.001  # for approximating DiracDelta distribution below

# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ("e", "eb", "x", "xb")


class Channel(BaseChannel):
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        # Note: `self.flavor` is set during __init__
        nu_flv = {'e': 12, 'eb': -12, 'x': 14, 'xb': -14}[self.flavor]

        evt = Event(2006012 if nu_flv > 0 else -2006012)
        evt.incoming_particles.append([nu_flv, eNu, 0, 0, 1])  # incoming nu
        evt.incoming_particles.append((6012, 11178, 0, 0, 1))  # carbon nucleus at rest
        evt.outgoing_particles.append([22, e_thr, dirx, diry, dirz])  # emitted gamma
        # evt.outgoing_particles.append([nu_flv, eNu-e_thr, 0, 0, 1])  # outgoing nu
        return evt

    # List with minimum & maximum energy of incoming neutrino.
    bounds_eNu = (e_thr, 100)

    def bounds_eE(self, eNu, *args):
        """Return kinematic bounds for integration over eE.

        Input:
            eNu:  neutrino energy (in MeV)
            args: [ignore this]
        Output:
            list with minimum & maximum allowed energy of outgoing (detected) particle
        """
        return [self.get_eE(eNu) - epsilon, self.get_eE(eNu) + epsilon]

    def get_eE(self, eNu, cosT=0):
        """Return energy (in MeV) of outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (in MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        return e_thr

    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        if eNu < e_thr or abs(self.get_eE(eNu) - eE) > epsilon:
            # This should never happen, since we set bounds for eE and eNu accordingly above
            # ... but just in case:
            return 0

        sigma = 1.08E-38 * ((eNu - e_thr) / 939)**2 * 1.11**2  # Assume beta * kappa = 1.11 (see discussion above)
        sigma *= (5.067731E10)**2  # convert cm^2 to MeV^-2: http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
        return sigma / (2 * epsilon)  # Ensure that integration over eE yields sigma

    def dSigma_dCosT(self, eNu, cosT):
        """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        # Energy dependence is unclear, so we use a constant value for now.
        if abs(cosT) > 1:
            return 0
        return 0.5

    def _bounds_eNu(self, eE):
        """Min/max neutrino energy that can produce a given positron energy."""
        return self.bounds_eNu
