"""Implementation of nu + 16O -> nu' + X. Where X is HERE NEED TO WRITE MORE HERE!!!!!

HERE NEED TO WRITE MORE HERE!!!!! - Talk about papers used, why and does it agree with T2K data?????

To determine the the differential cross section dSigma/dE (eNu, eE) from the
total cross section, we approximate a DiracDelta function with one that is
2*epsilon wide and 1/(2*epsilon) high, so that the integral is 1.
"""

from sntools.event import Event
from sntools.interaction_channels import BaseChannel
from scipy.interpolate import interp1d
import numpy as np
import random

""" The following table of data is taken from Suzuki et al. 2018. 
This is included to generate a spline function of the cross-section as a function
of the energy of the incoming neutrino
"""

x_sec_data = [(14  , None    , None    , 0.0     , None    , None    , None    , None    , None    , None    , None    , 0.0     ),
              (15.0, 1.35E-07, None    , 6.05E-04, None    , None    , None    , None    , 1.35E-25, None    , None    , 6.05E-04), 
              (16.0, 5E-07   , None    , None    , None    , None    , None    , None    , 1E-19   , None    , None    , None    ), #Estimated
              (16.5, None    , None    , 3E-03   , None    , None    , None    , None    , None    , None    , None    , None    ), #Estimated
              (18.0, None    , None    , None    , None    , None    , None    , None    , None    , None    , None    , 1.00E-02), #Estimated
              (18.5, None    , 0.0     , None    , None    , None    , None    , None    , None    , None    , None    , None    ), #Estimated
              (19.5, None    , 1E-04   , None    , None    , None    , None    , None    , None    , None    , None    , None    ), #Estimated
              (20.0, 1.02E-05, 1.93E-04, 3.08E-02, 0.0     , 0.0     , 0.0     , 0.0     , 1.47E-10, None    , 0.0     , 3.10E-02), 
              (25.0, 9.97E-05, 1.63E-02, 1.98E-01, 2.88E-05, 1.09E-05, 2.46E-09, 1.82E-05, 1.74E-06, 0.0     , 4.30E-08, 2.14E-01), 
              (26.0, None    , None    , None    , None    , 1E-04   , 1E-06   , None    , None    , None    , None    , None    ), #Estimated
              (27.0, None    , None    , None    , 6E-04   , None    , None    , 4E-04   , None    , None    , 1E-04   , None    ), #Estimated
              (28.8, None    , None    , None    , None    , None    , None    , None    , None    , 1E-05   , None    , None    ), #Estimated
              (30.0, 5.54E-04, 1.29E-01, 7.76E-01, 7.92E-03, 5.71E-03, 1.82E-03, 4.19E-03, 4.14E-03, 3.76E-05, 4.83E-03, 9.34E-01), 
              (35.0, 1.83E-03, 4.80E-01, 2.16E+00, 6.43E-02, 4.74E-02, 1.48E-02, 3.56E-02, 3.17E-02, 6.86E-04, 5.39E-02, 2.89E+00), 
              (40.0, 4.46E-03, 1.22E+00, 4.76E+00, 2.06E-01, 1.51E-01, 4.52E-02, 1.18E-01, 9.55E-02, 3.74E-03, 1.99E-01, 6.80E+00), 
              (45.0, 9.00E-03, 2.49E+00, 9.01E+00, 4.68E-01, 3.41E-01, 9.96E-02, 2.75E-01, 2.07E-01, 1.09E-02, 4.96E-01, 1.34E+01), 
              (50.0, 1.60E-02, 4.44E+00, 1.52E+01, 8.87E-01, 6.44E-01, 1.85E-01, 5.36E-01, 3.81E-01, 2.43E-02, 9.97E-01, 2.33E+01), 
              (55.0, 2.60E-02, 7.17E+00, 2.37E+01, 1.49E+00, 1.08E+00, 3.06E-01, 9.23E-01, 6.27E-01, 4.53E-02, 1.75E+00, 3.71E+01), 
              (60.0, 3.94E-02, 1.07E+01, 3.45E+01, 2.31E+00, 1.68E+00, 4.69E-01, 1.45E+00, 9.57E-01, 7.48E-02, 2.80E+00, 5.50E+01), 
              (65.0, 5.66E-02, 1.52E+01, 4.77E+01, 3.36E+00, 2.43E+00, 6.75E-01, 2.14E+00, 1.37E+00, 1.13E-01, 4.17E+00, 7.72E+01), 
              (70.0, 7.78E-02, 2.04E+01, 6.29E+01, 6.64E+00, 3.35E+00, 9.24E-01, 2.97E+00, 1.87E+00, 1.60E-01, 5.85E+00, 1.05E+02), 
              (80.0, 1.34E-01, 3.28E+01, 9.86E+01, 7.82E+00, 5.63E+00, 1.54E+00, 5.06E+00, 3.11E+00, 2.78E-01, 1.01E+01, 1.65E+02), 
              (90.0, 2.10E-01, 4.68E+01, 1.38E+02, 1.17E+01, 8.36E+00, 2.27E+00, 7.57E+00, 4.58E+00, 4.23E-01, 1.53E+01, 2.35E+02)] 

gam, n, p, d, pp, H3, He3, alp, alp_n, alp_p, tot=[],[],[],[],[],[],[],[],[],[],[] 
a = (gam, n, p, d, pp, H3, He3, alp, alp_n, alp_p, tot)
for j in range(1, 12):   # sorting through the data removing any 'None' data points and returning
    e_Nu = []
    s = []     # a list of two lists  of the form [energy, cross-section for given final state]
    for i in range(0, len(x_sec_data)):
        if x_sec_data[i][j] is None:    continue
        else:
            s.append(x_sec_data[i][j])
            e_Nu.append(x_sec_data[i][0])
    a[j-1].append(e_Nu)
    a[j-1].append(s)
    
data = {'gamma':gam, 'neutron':n, 'proton':p, 'deuteron':d,     # Dictionary of data sets for the final states
        'pp':pp, 'H3':H3, 'He3':He3, 'alpha':alp, 
        'alpha_n':alp_n, 'alpha_p':alp_p, 'total':tot}

def fit(state):
    eNu = data[state][0]
    f = interp1d(eNu, data[state][1], kind='cubic', fill_value='extrapolate', bounds_error = False)
    return f

e_thr = 18.7  # approximate energy threshold of neutron knock-out reaction with no gamma emission (MeV)
e_thr_g = e_thr + 6.18  # energy threshold of neutron knock-out reaction with gamma emission (MeV)
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

        eE = self.get_eE(eNu, dirz)

        evt = Event(2006012 if nu_flv > 0 else -2006012)
        evt.incoming_particles.append([nu_flv, eNu, 0, 0, 1])  # incoming nu
        evt.incoming_particles.append((8016, 14900, 0, 0, 1))  # oxygen-16 nucleus at rest

        if self.gamma(eNu) is True:   # excited O-16 nucleus emits a neutron and a photon
            g_dirx, g_diry, g_dirzN = self.get_photon_direction()
            evt.outgoing_particles.append([22, 6.18, g_dirx, g_diry, g_dirzN])  # emitted photon
            evt.outgoing_particles.append([2112, eE, dirx, diry, dirz])  # emitted neutron
        else:       # excited O-16 nucleus emits just a neutron
            evt.outgoing_particles.append([2112, eE, dirx, diry, dirz])  # emitted neutron
        """
        NOTE: Need to define gamma and neutron energies as de-excitation energy and distribution energy respectively
        """
        # evt.outgoing_particles.append([nu_flv, eNu-e_thr, 0, 0, 1])  # outgoing nu
        return evt

    def gamma(self, eNu):
        """Returns True if a gamma emission occurs.
        If eNu > e_thr_g, then the O-15 nucleus is left in the excited state in ~70% of interactions.
        Therefore a de-excitation gamma is produced in addition to the knocked-out neutron.
        """
        r = random.random()
        if eNu > e_thr_g and r > 0.310107949:
            return True
        else: 
            return False

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
        if self.gamma() is True:
            eE = random.random()*(eNu - e_thr_g) # energy of emitted neutron generated randomly from energy excess of neutrino over gamma emission threshold energy
        else: eE = random.random()*(eNu-e_thr)    # energy of emitted neutron generated randomly from energy excess of neutrino over threshold energy           
        return eE
    
    def get_photon_direction(self):
        """
        In the case where nucleon emission leaves the daughter nucleus in an
        excited state, a photon is emitted to de-excite the nucleus.
        
        This function returns a random direction (normalised to 1) for this 
        photon emission.
        """
        phi = 2*np.pi*random.random()
        theta = np.pi*random.random()
        
        dirx = np.sin(theta)*np.cos(phi)
        diry = np.sin(theta)*np.sin(phi)
        dirz = np.cos(theta)
        
        return (dirx,diry,dirz)
    
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

        sigma = fit('neutron')(eNu)  # cross-section at eNu from the interp1d fit of Suzuki et al. 2018 data
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
    
    # List with minimum & maximum energy of incoming neutrino.
    bounds_eNu = (e_thr, 100)

    def _bounds_eNu(self, eE):
        """Min/max neutrino energy that can produce a given positron energy."""
        return self.bounds_eNu
