"""Implementation of:
    nu + 16O -> nu' + 16O*, 16O* -> 15O + n.

Based on data provided in Suzuki et al. 2018 (Phys. Rev. C 98, 034613).
A spline fit is used to obtain cross-section as a function of neutrino energy.
The threshold energy of the interaction is estimated from a by-eye fit as no 
explicit threshold energy could be found.

To determine the the differential cross section dSigma/dE (eNu, eE) from the
total cross section, we approximate a DiracDelta function with one that is
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""

from sntools.event import Event
from sntools.interaction_channels import BaseChannel
from scipy.interpolate import interp1d
import random

e_thr = 19.5  # approximate energy threshold of neutron emission
mN = 939.7 # neutron mass (MeV)
epsilon = 0.001  # for approximating DiracDelta distribution below

# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ("e", "eb", "x", "xb")

data = [[e_thr, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 80.0, 90.0],
 [0.0, 0.000193, 0.0163, 0.129, 0.48, 1.22, 2.49, 4.44, 7.17, 10.7, 15.2, 20.4, 32.8, 46.8]]
 #list of energies and cross-sections for o16nc neutron emission taken from Suzuki et al 2018.

# Spline fit of the partial cross-section as a function of energy
fit = interp1d(data[0], data[1], kind='cubic', fill_value='extrapolate', bounds_error = False)

class Channel(BaseChannel):
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        # Note: `self.flavor` is set during __init__
        nu_flv = {'e': 12, 'eb': -12, 'x': 14, 'xb': -14}[self.flavor]

        eE = self.get_eE(eNu, dirz)

        evt = Event(2008016 if nu_flv > 0 else -2008016)
        evt.incoming_particles.append([nu_flv, eNu, 0, 0, 1])  # incoming nu
        evt.incoming_particles.append((8016, 14900, 0, 0, 1))  # oxygen-16 nucleus at rest
        evt.outgoing_particles.append([2112, eE, dirx, diry, dirz])  # emitted neutron
        # evt.outgoing_particles.append([nu_flv, eNu-e_thr, 0, 0, 1])  # outgoing nu
        return evt
    
    def bounds_eE(self, eNu, *args):
        """Return kinematic bounds for integration over eE.

        Input:
            eNu:  neutrino energy (in MeV)
            args: [ignore this]
        Output:
            list with minimum & maximum allowed energy of outgoing (detected) particle
        """
        eE = self.get_eE(eNu)
        return [eE - epsilon, eE + epsilon]
        
    def get_eE(self, eNu, cosT=0):
        """Return energy (in MeV) of outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (in MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        eE = random.random()*(eNu-e_thr) + mN       
        return eE
    
    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        if eNu < e_thr:
            # This should never happen, since we set bounds for eNu accordingly above
            # ... but just in case:
            return 0

        sigma = fit(eNu)*10**(-42)  # cross-section at eNu from the interp1d fit 
                                                # of Suzuki et al. 2018 data (units cm^2)
        sigma *= (5.067731E10)**2  # convert cm^2 to MeV^-2: 
                   # http://www.wolframalpha.com/input/?i=cm%2F(hbar+*+c)+in+MeV%5E(-1)
        return sigma / (2 * epsilon)  # Ensure that integration over eE yields sigma

    def dSigma_dCosT(self, eNu, cosT):
        """Return differential cross section in MeV^-2 as a function of the 
        emission angle of the outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        # Energy dependence is unclear, so we use a constant value for now.
        if abs(cosT) > 1:
            return 0
        return 0.5
    
    # List with minimum & maximum energy of incoming neutrino.
    bounds_eNu = (e_thr, 100)

    def _bounds_eNu(self, eE):
        """Min/max neutrino energy that can produce a given neutron energy."""
        return (e_thr + eE - mN, 100)