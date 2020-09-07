"""Implementation of anti-nu_e + p -> n + e+

Based on Strumia/Vissani (2003), arXiv:astro-ph/0302055.
"""

from __future__ import division

from math import pi, sqrt, log

from sntools.event import Event


def generate_event(eNu, dirx, diry, dirz):
    """Return an event with the appropriate incoming/outgoing particles.

    Input:
        eNu: neutrino energy
        dirx, diry, dirz: direction of outgoing particle (normalized to 1)
    """
    eE = get_eE(eNu, dirz)
    eN, dirxN, diryN, dirzN = get_neutron_kinematics(eNu, eE, dirx, diry, dirz)

    evt = Event(-1001001)
    evt.incoming_particles.append((-12, eNu, 0, 0, 1))  # incoming neutrino
    evt.incoming_particles.append((2212, mP, 0, 0, 1))  # proton at rest
    evt.outgoing_particles.append((-11, eE, dirx, diry, dirz))  # outgoing positron
    evt.outgoing_particles.append((2112, eN, dirxN, diryN, dirzN))  # outgoing neutron
    return evt


# List of neutrino flavors ("e", "eb", "x", "xb") that interact in this channel.
possible_flavors = ["eb"]

mN = 939.56563  # neutron mass (MeV)
mP = 938.27231  # proton mass (MeV)
mE = 0.5109907  # electron mass (MeV)
mPi = 139.56995  # pion mass (MeV)
delta = mN - mP
mAvg = (mP + mN) / 2
alpha = 1 / 137.035989  # fine structure constant
gF = 1.16639e-11  # Fermi coupling constant
sigma0 = 2 * mP * gF**2 * 0.9746**2 / (8 * pi * mP**2)  # from eqs. (3), (11)


def dSigma_dE(eNu, eE):  # eqs. (11), (3)
    """Return differential cross section in MeV^-2.

    Inputs:
        eNu: neutrino energy
        eE:  energy of outgoing (detected) particle
    """
    if eNu < eThr or eE < bounds_eE(eNu)[0] or eE > bounds_eE(eNu)[1]:
        return 0
    # above eq. (11)
    s_minus_u = 2 * mP * (eNu + eE) - mE**2
    t = mN**2 - mP**2 - 2 * mP * (eNu - eE)

    # eq. (7)
    x = 0 + t / (4 * mAvg**2)
    y = 1 - t / 710**2
    z = 1 - t / 1030**2
    f1 = (1 - 4.706 * x) / ((1 - x) * y**2)
    f2 = 3.706 / ((1 - x) * y**2)
    g1 = -1.27 / z**2
    g2 = 2 * g1 * mAvg**2 / (mPi**2 - t)

    A = (t - mE**2) * (
            4 * f1**2 * (4 * mAvg**2 + t + mE**2)
            + 4 * g1**2 * (-4 * mAvg**2 + t + mE**2)
            + f2**2 * (t**2 / mAvg**2 + 4 * t + 4 * mE**2)
            + 4 * mE**2 * t * g2**2 / mAvg**2
            + 8 * f1 * f2 * (2 * t + mE**2)
            + 16 * mE**2 * g1 * g2) \
        - delta**2 * (
            (4 * f1**2 + t * f2**2 / mAvg**2) * (4 * mAvg**2 + t - mE**2)
            + 4 * g1**2 * (4 * mAvg**2 - t + mE**2)
            + 4 * mE**2 * g2**2 * (t - mE**2) / mAvg**2
            + 8 * f1 * f2 * (2 * t - mE**2)
            + 16 * mE**2 * g1 * g2) \
        - 32 * mE**2 * mAvg * delta * g1 * (f1 + f2)
    A /= 16

    B = 16 * t * g1 * (f1 + f2) + 4 * mE**2 * delta * (f2**2 + f1 * f2 + 2 * g1 * g2) / mAvg
    B /= 16

    C = 4 * (f1**2 + g1**2) - t * f2**2 / mAvg**2
    C /= 16

    abs_M_squared = A - B * s_minus_u + C * s_minus_u**2  # eq. (5)
    rad_correction = alpha / pi * (6.00352 + 3 / 2 * log(mP / (2 * eE)) + 1.2 * (mE / eE)**1.5)  # eq. (14)

    result = sigma0 / eNu**2 * abs_M_squared * (1 + rad_correction)

    if result < 0:
        raise ValueError("Calculated negative cross section for E_nu=%f, E_e=%f. Aborting..." % (eNu, eE))

    return result


def dSigma_dCosT(eNu, cosT):  # eq. (20)
    """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    epsilon = eNu / mP
    eE = get_eE(eNu, cosT)
    pE = sqrt(eE**2 - mE**2)
    dE_dCosT = pE * epsilon / (1 + epsilon * (1 - cosT * eE / pE))
    return dE_dCosT * dSigma_dE(eNu, eE)


def get_eE(eNu, cosT):  # eq. (21)
    """Return energy (in MeV) of outgoing (detected) particle.

    Input:
        eNu:  neutrino energy (in MeV)
        cosT: cosine of the angle between neutrino and outgoing (detected) particle
    """
    epsilon = eNu / mP
    kappa = (1 + epsilon)**2 - (epsilon * cosT)**2
    return ((eNu - delta_cm) * (1 + epsilon) + epsilon * cosT * sqrt((eNu - delta_cm)**2 - mE**2 * kappa)) / kappa


def get_neutron_kinematics(eNu, eE, dirx, diry, dirz):
    eN = mP + eNu - eE  # neutron energy

    # calculate 3-momentum of neutron ...
    pE = sqrt(eE**2 - mE**2)
    pN_x = - dirx * pE
    pN_y = - diry * pE
    pN_z = eNu - dirz * pE
    # ... and normalize it to get unit vector of neutron direction
    pN = sqrt(pN_x**2 + pN_y**2 + pN_z**2)
    return (eN, pN_x / pN, pN_y / pN, pN_z / pN)


# Bounds for integration over eE
delta_cm = (mN**2 - mP**2 - mE**2) / (2 * mP)


def bounds_eE(eNu, *args):
    """Return kinematic bounds for integration over eE.

    Input:
        eNu:  neutrino energy (in MeV)
        args: [ignore this]
    Output:
        list with minimum & maximum allowed energy of outgoing (detected) particle
    """
    s = 2 * mP * eNu + mP**2
    pE_cm = sqrt((s - (mN - mE)**2) * (s - (mN + mE)**2)) / (2 * sqrt(s))
    eE_cm = (s - mN**2 + mE**2) / (2 * sqrt(s))

    eE_min = eNu - delta_cm - eNu / sqrt(s) * (eE_cm + pE_cm)
    eE_max = eNu - delta_cm - eNu / sqrt(s) * (eE_cm - pE_cm)
    return [eE_min, eE_max]


# Bounds for integration over eNu
eThr = ((mN + mE)**2 - mP**2) / (2 * mP)  # threshold energy for IBD: ca. 1.8 MeV
bounds_eNu = [eThr, 100]


def _bounds_eNu(eE):
    """Min/max neutrino energy that can produce a given positron energy.

    This is an approximation to simplify numerical integration; the precise range has to be enforced separately.
    """
    eNu_min = eE + delta_cm
    eNu_max = eNu_min / (1 - 2 * eNu_min / mP)
    return (eNu_min, eNu_max)
