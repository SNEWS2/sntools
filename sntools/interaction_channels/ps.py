"""
Implementation of nu + p -> nu + p

Largely based on John F Beacom et al., 2002: 10.1103/PhysRevD.66.033001.
This neutrino-proton scattering process only has a NC channel. The calculation assumes neutrinos at SN neutrino energies.

Written by S. Valder (2022), based on other interaction channels in sntools.
"""

from math import pi, sqrt, log
from scipy import integrate

from sntools.event import Event
from sntools.interaction_channels import BaseChannel, cherenkov_threshold


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["e", "eb", "x", "xb"]

sin2theta_w = 0.23155  # weak mixing angle
mP = 938.27205  # proton mass (MeV)
gF = 1.16637e-11  # Fermi coupling constant
gA_0 = 1.267 # axial proton form factor; from Beringer et al. 2012, doi: 10.1103/PhysRevD.86.010001
eta = 0.12 # proton strangeness; from Ahrens et al. 1987, doi 10.1103/PhysRevD.35.785
cA_0 = (gA_0 * (1 + eta)) / 2 # axial vector coupling constant
cV = (1 - (4 * sin2theta_w)) / 2 # axial vector coupling constant


class Channel(BaseChannel):
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        incoming_flv = {'e': 12, 'eb': -12, 'x': 14, 'xb': -14}[self.flavor]
        eE = self.get_eE(eNu, dirz)

        evt = Event(2001001 if incoming_flv > 0 else -2001001)
        evt.incoming_particles.append((incoming_flv, eNu, 0, 0, 1))  # incoming neutrino
        evt.incoming_particles.append((2212, mP, 0, 0, 1))  # proton at rest
        evt.outgoing_particles.append((2212, eE, dirx, diry, dirz))  # outgoing proton
        return evt
    
    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        if eE < self.bounds_eE(eNu)[0] or eE > self.bounds_eE(eNu)[1]:
            return 0
 
        if self.flavor in ("e", "x"):  # neutrino
            cA = cA_0

        if self.flavor in ("eb", "xb"):  # antineutrino
            cA = -cA_0

        T = eE - mP # Kinetic energy of outgoing proton

        result = (((gF**2 * mP) / (2 * pi * eNu**2)) * (((cV + cA)**2 * eNu**2) + ((cV - cA)**2 * (eNu - T)**2) - ((cV**2 - cA**2) * mP * T)))

        if result < 0:
            raise ValueError(f"Calculated negative cross section for E_nu={eNu}, E_e={eE}. Aborting...")
        
        return result

    def get_eE(self, eNu, cosT):
        """Return energy of scattered proton (in MeV).

        Input:
            eNu:  neutrino energy (in MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """

        return (((2 * cosT**2 * eNu**2) / (mP)) + mP)  # calculation from 10.1103/PhysRevD.66.033001, first order approximation


    def dSigma_dCosT(self, eNu, cosT):
        """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        if cosT < 0:  # backward scattering is kinematically impossible
            return 0

        dE_dCosT = (4 * cosT * eNu**2) / mP # differentiate eE wrt cosT
        eE = self.get_eE(eNu, cosT)
        return dE_dCosT * self.dSigma_dE(eNu, eE) 

    # Minuimum energy of the outgoing proton is the mass of the proton
    eE_min = mP

    def bounds_eE(self, eNu, *args):  # ignore additional arguments handed over by integrate.nquad()
        """Return kinematic bounds for integration over eE.

        Input:
            eNu:  neutrino energy (in MeV)
            args: [ignore this]
        Output:
            list with minimum & maximum allowed energy of outgoing (detected) particle
        """
        
        #eE_max = mP + ((2 * eNu**2) / mP)   # this is get_eE(eNu, cosT=1); could use 
        eE_max = mP + ((2 * eNu**2) / (mP + 2 * eNu))

        return [self.eE_min, eE_max]

    # Bounds for integration over eNu
    def eNu_min(self, eE):
        T = eE - mP  
        return ((T + (sqrt(T * (T + (2 * mP))))) / 2) # see 10.1103/PhysRevD.66.033001

    bounds_eNu = [eNu_min(None, eE_min), 100]

    def _bounds_eNu(self, eE):
        """Min/max neutrino energy that can produce a given positron energy."""
        return [self.eNu_min(eE), self.bounds_eNu[1]]

