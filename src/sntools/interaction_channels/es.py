"""Implementation of nu + e- -> nu + e-

Based on the Appendices of Bahcall et al. 1995 (https://doi.org/10.1103/PhysRevD.51.6146).
This calculation includes radiative corrections from QCD & QED effects.

For differences between neutrinos/antineutrinos and a derivation of directionality,
see https://www.kvi.nl/~loehner/saf_seminar/2010/neutrino-electron-interactions.pdf
Note that it uses different conventions (e.g. minus signs) from Bahcall et al.!
"""

from math import pi, sqrt, log
from scipy import integrate

from sntools.event import Event
from sntools.interaction_channels import BaseChannel, cherenkov_threshold


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["e", "eb", "x", "xb"]

sin2theta_w = 0.2317  # weak mixing angle
alpha = 1 / 137.036  # fine structure constant
mE = 0.5109989  # electron mass (MeV)
gF = 1.16637e-11  # Fermi coupling constant
rho_NC = 1.0126  # numerical factor from Bahcall et al.

_cache = {}  # save time by avoiding repeat calculations
def spence(n):
    n = round(n, 6)  # negligible decrease in accuracy to save time & memory
    if n not in _cache:
        _cache[n] = integrate.quad(lambda t: log(abs(1 - t)) / t, 0, n)[0]
    return _cache[n]


class Channel(BaseChannel):
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        incoming_flv = {'e': 12, 'eb': -12, 'x': 14, 'xb': -14}[self.flavor]
        eE = self.get_eE(eNu, dirz)

        evt = Event(98 if incoming_flv > 0 else -98)
        evt.incoming_particles.append((incoming_flv, eNu, 0, 0, 1))  # incoming neutrino
        evt.incoming_particles.append((11, mE, 0, 0, 1))  # electron at rest
        evt.outgoing_particles.append((11, eE, dirx, diry, dirz))  # outgoing electron
        return evt

    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        if eE < self.bounds_eE(eNu)[0] or eE > self.bounds_eE(eNu)[1]:
            return 0

        # Appendix A: Radiative Corrections
        L = sqrt(eE**2 - mE**2)
        beta = L / eE
        T = eE - mE  # kinetic energy of recoil electron
        z = T / eNu
        x = sqrt(1 + 2 * mE / T)
        i = 1 / 6 * (1 / 3 + (3 - x**2) * (x / 2 * log((x + 1) / (x - 1)) - 1))

        if self.flavor in ("e", "eb"):
            k = 0.9791 + 0.0097 * i
        elif self.flavor in ("x", "xb"):
            k = 0.9970 - 0.00037 * i

        g1 = rho_NC * (0.5 - k * sin2theta_w)
        g2 = -rho_NC * k * sin2theta_w

        if self.flavor == "e":
            gL = g1 - 1
            gR = g2
        elif self.flavor == "eb":
            gL = g2
            gR = g1 - 1
        elif self.flavor == "x":
            gL = g1
            gR = g2
        elif self.flavor == "xb":
            gL = g2
            gR = g1

        # Appendix B: QED Effects
        f0 = eE / L * log((eE + L) / mE) - 1  # common factor of all three f_*
        log_zmE = log(1 - z - mE / (eE + L))  # Warning: imprecise at low E, throws ValueError in extreme cases
        log1z = log(1 - z)

        # fMinus(z)
        f1 = f0 * (2 * log_zmE - log1z - log(z) / 2 - 5 / 12) \
            + 0.5 * (spence(z) - spence(beta)) \
            - 0.5 * log1z**2 - (11 / 12 + z / 2) * log1z \
            + z * (log(z) + 0.5 * log(2 * eNu / mE)) \
            - (31 / 18 + 1 / 12 * log(z)) * beta \
            - 11 / 12 * z + z**2 / 24

        # (1 - z)**2 * fPlus(z)
        f2 = f0 * ((1 - z)**2 * (2 * log_zmE - log1z - log(z) / 2 - 2 / 3) - (z**2 * log(z) + 1 - z) / 2) \
            - (1 - z)**2 / 2 * (log1z**2 + beta * (spence(1 - z) - log(z) * log1z)) \
            + log1z * (z**2 / 2 * log(z) + (1 - z) / 3 * (2 * z - 0.5)) \
            - z**2 / 2 * spence(1 - z) - z * (1 - 2 * z) / 3 * log(z) - z * (1 - z) / 6 \
            - beta / 12 * (log(z) + (1 - z) * (115 - 109 * z) / 6)

        # fPlusMinus(z)
        f3 = f0 * 2 * log_zmE

        result = 2 * mE * gF**2 / pi * (gL**2 * (1 + alpha / pi * f1)
                                        + gR**2 * ((1 - z)**2 + alpha / pi * f2)
                                        - gR * gL * mE / eNu * z * (1 + alpha / pi * f3))
        if result < 0:
            if eNu < cherenkov_threshold:
                # Approximations in f_* may be imprecise at very low energies.
                # This is below threshold in HK anyway, so we suppress it.
                result = 0
            else:
                raise ValueError(f"Calculated negative cross section for E_nu={eNu}, E_e={eE}. Aborting...")

        return result

    def get_eE(self, eNu, cosT):
        """Return energy of scattered electron (in MeV).

        Input:
            eNu:  neutrino energy (in MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        return mE + (2 * mE * eNu**2 * cosT**2) / ((mE + eNu)**2 - eNu**2 * cosT**2)

    def dSigma_dCosT(self, eNu, cosT):
        """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        if cosT < 0:  # backward scattering is kinematically impossible
            return 0

        dE_dCosT = 4 * mE * eNu**2 * (mE + eNu)**2 * cosT / ((mE + eNu)**2 - eNu**2 * cosT**2)**2
        eE = self.get_eE(eNu, cosT)
        return dE_dCosT * self.dSigma_dE(eNu, eE)

    eE_min = cherenkov_threshold

    def bounds_eE(self, eNu, *args):  # ignore additional arguments handed over by integrate.nquad()
        """Return kinematic bounds for integration over eE.

        Input:
            eNu:  neutrino energy (in MeV)
            args: [ignore this]
        Output:
            list with minimum & maximum allowed energy of outgoing (detected) particle
        """
        eE_max = mE + 2 * eNu**2 / (2 * eNu + mE)  # this is get_eE(eNu, cosT=1)
        return [self.eE_min, eE_max]

    # Bounds for integration over eNu
    def eNu_min(self, eE):
        T = eE - mE
        return T / 2. * (1 + sqrt(1 + 2 * mE / T))  # inversion of eE_max(eNu)

    bounds_eNu = [eNu_min(None, eE_min), 100]

    def _bounds_eNu(self, eE):
        """Min/max neutrino energy that can produce a given positron energy."""
        return [self.eNu_min(eE), self.bounds_eNu[1]]
