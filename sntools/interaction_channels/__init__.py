from abc import ABC, abstractmethod


cherenkov_threshold = 0.77  # Cherenkov threshold of electron/positron in water


class BaseChannel(ABC):
    """Abstract base class that defines the interface for all submodules in `sntools.interaction_channels`."""
    def __init__(self, flv) -> None:
        self.flavor = flv

    @abstractmethod
    def generate_event(self, eNu, dirx, diry, dirz):
        """Return an event with the appropriate incoming/outgoing particles.

        Input:
            eNu: neutrino energy
            dirx, diry, dirz: direction of outgoing particle (normalized to 1)
        """
        pass

    @abstractmethod
    def dSigma_dE(self, eNu, eE):
        """Return differential cross section in MeV^-2.

        Inputs:
            eNu: neutrino energy
            eE:  energy of outgoing (detected) particle
        """
        pass

    @abstractmethod
    def dSigma_dCosT(eNu, cosT):
        """Return differential cross section in MeV^-2 as a function of the emission angle of the outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        pass

    @abstractmethod
    def get_eE(self, eNu, cosT):
        """Return energy (in MeV) of outgoing (detected) particle.

        Input:
            eNu:  neutrino energy (in MeV)
            cosT: cosine of the angle between neutrino and outgoing (detected) particle
        """
        pass

    @abstractmethod
    def bounds_eE(self, eNu, *args):
        """Return kinematic bounds for integration over eE.

        Input:
            eNu:  neutrino energy (in MeV)
            args: [ignore this]
        Output:
            list with minimum & maximum allowed energy of outgoing (detected) particle
        """
        pass

    @property
    @abstractmethod
    def bounds_eNu(self):
        """Tuple with minimum & maximum energy of incoming neutrino."""
        return (0, 100)
