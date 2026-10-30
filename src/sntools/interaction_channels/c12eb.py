"""Implementation of anti-nu_e + 12C -> X + e+

Based on arXiv:nucl-th/9903022, Table V. Following arXiv:1204.4231 (eqn. III.1),
we fit a power series to the tabulated values. Furthermore, we ignore excited
states of the final nucleus, so that the observed energy is eE = eNu - e_thr.

To determine the the differential cross section dSigma/dE (eNu, eE) from the
total cross section, we approximate a DiracDelta function with one that is
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""

from sntools.event import Event
from sntools.interaction_channels import BaseChannel, cherenkov_threshold

import numpy as np
from scipy.interpolate import PchipInterpolator

e_thr = 14.39  # energy threshold of this reaction (arXiv:1507.05613, p. 74)
epsilon = 0.001  # for approximating DiracDelta distribution below

#Table V from arXiv:nucl-th/9903022, plus extra entry for e_thr
#E_nu in MeV, sigma in cm^2.
_Enu_arr_interp = np.array([
    e_thr,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,  80.0,  90.0, 100.0,
    110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
    210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0])
_sigma_arr_interp = 1e-42*np.array([
    0.000,   0.8,  6.05,  18.4,  41.4,  78.0, 130.0, 199.0, 284.0, 386.0,
    501.0, 629.0, 768.0, 917.0, 1080., 1240., 1420., 1590., 1780., 1960.,
    2150., 2340., 2530., 2720., 2910., 3110., 3300., 3500., 3700., 3900.])

_sigma_interp = PchipInterpolator(_Enu_arr_interp, _sigma_arr_interp, extrapolate=False)

# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["eb"]


class Channel(BaseChannel):
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        eE = self.get_eE(eNu, dirz)
        evt = Event(-1006012)
        evt.incoming_particles.append([-12, eNu, 0, 0, 1])  # incoming neutrino
        evt.incoming_particles.append((6012, 11178, 0, 0, 1))  # carbon nucleus at rest
        evt.outgoing_particles.append([-11, eE, dirx, diry, dirz])  # outgoing positron
        return evt

    # List with minimum & maximum energy of incoming neutrino.
    bounds_eNu = [e_thr + cherenkov_threshold, 100]

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
        return eNu - e_thr

    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        if abs(self.get_eE(eNu) - eE) > epsilon:
            # This should never happen, since we set bounds_eE() accordingly above
            # ... but just in case:
            return 0

        #Polynomial fit: Do not use, since it becomes negative near e_thr.
        #a1, a2, a3 = -106.3, 25.15, 0.3697
        #sigma = 9.55e-46 * (a1 * (eNu - e_thr) + a2 * (eNu - e_thr)**2 + a3 * (eNu - e_thr)**3)
        
        #Ensure zero cross section below threshold.
        if eNu < e_thr:
            return 0
        
        sigma = _sigma_interp(eNu)
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
        return (eE + e_thr - epsilon, eE + e_thr + epsilon)
