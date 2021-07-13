from abc import ABC, abstractmethod
from importlib import import_module
from math import ceil, floor


class BaseFlux(ABC):
    """Abstract base class that defines the interface for all submodules in `sntools.formats`."""
    @abstractmethod
    def parse_input():
        pass

    @abstractmethod
    def prepare_evt_gen():
        pass

    @abstractmethod
    def nu_emission():
        pass


class WeightedFlux(BaseFlux):
    def __init__(self, unweighted_flux, weight) -> None:
        super().__init__()
        self._flux = unweighted_flux
        self._weight = weight

    def __repr__(self) -> str:
        return f"WeightedFlux({self._flux}, {self._weight})"

    def __getattr__(self, name: str):
        if name != "_flux":
            return getattr(self._flux, name)

    def parse_input(self, *args):
        self._flux.parse_input(*args)

    def prepare_evt_gen(self, *args):
        self._flux.prepare_evt_gen(*args)

    def nu_emission(self, eNu, time):
        return self._weight * self._flux.nu_emission(eNu, time)


class CompositeFlux:
    """Neutrino flux encompassing different BaseFlux components; zero, one or more for each neutrino flavor."""

    def __init__(self) -> None:
        # Dictionary with up to 4 components (e, eb, x, xb), each of which is a list of BaseFlux instances.
        self.components = {'e': [], 'eb': [], 'x': [], 'xb': []}
        self._repr = super().__repr__()

    def __repr__(self) -> str:
        return self._repr

    @classmethod
    def from_file(cls, file, format, starttime=None, endtime=None):
        """Create a CompositeFlux from an input file."""
        self = CompositeFlux()
        self._repr = f"CompositeFlux.from_file('{file}', '{format}', {starttime}, {endtime})"
        format = import_module('sntools.formats.' + format)

        for flv in ('e', 'eb', 'x', 'xb'):
            f = format.Flux()
            f.parse_input(file, flv, starttime, endtime)
            self.components[flv] = [f]

        return self

    def transformed_by(self, transformation, distance):
        """Apply transformation and 1/d^2 scaling to calculate flux at detector."""
        tf = CompositeFlux()
        tf._repr = f"{self.__repr__()}.transformed_by(transformation={transformation}, distance={distance})"

        for detected_flv in tf.components:
            for (original_flv, scale) in transformation.components_producing(detected_flv):
                for flux in self.components[original_flv]:
                    wf = WeightedFlux(flux, scale * (10.0 / distance)**2)
                    tf.components[detected_flv].append(wf)

        return tf


class SNEWPYFlux(BaseFlux):
    """Adapter class to turn a SNEWPY.models.SupernovaModel component into an sntools.formats.BaseFlux"""

    def __init__(self, sn_model, flv, starttime, endtime) -> None:
        from astropy import units as u
        from snewpy.neutrino import Flavor

        super().__init__()
        self._flv = {'e': Flavor.NU_E, 'eb': Flavor.NU_E_BAR, 'x': Flavor.NU_X, 'xb': Flavor.NU_X_BAR}[flv]
        self._sn_model = sn_model

        times = [t.to(u.ms).value for t in sn_model.get_time()]
        self.starttime = get_starttime(starttime, times[0])
        self.endtime = get_endtime(endtime, times[-1])
        self.raw_times = times  # TODO: enforce starttime/endtime in self.raw_times

    def parse_input(self, *args):
        # handled in SNEWPYCompositeFlux.from_file()
        pass

    def prepare_evt_gen(self, *args):
        pass

    def nu_emission(self, eNu, time):
        from astropy import units as u
        nl = self._sn_model.get_initialspectra(time * u.ms, eNu * u.MeV)[self._flv]
        # SNEWPY uses cgs units internally, so convert before returning
        return nl.to(1 / u.MeV / u.ms).value


class SNEWPYCompositeFlux(CompositeFlux):
    """Adapter class to turn a SNEWPY.models.SupernovaModel into an sntools.formats.CompositeFlux"""

    @classmethod
    def from_file(cls, file, format, starttime=None, endtime=None):
        """Create a SNEWPYCompositeFlux from an input file."""
        self = SNEWPYCompositeFlux()
        self._repr = f"SNEWPYCompositeFlux.from_file('{file}', '{format}', {starttime}, {endtime})"

        sn_model = getattr(import_module('snewpy.models'), format)(file)

        for flv in ('e', 'eb', 'x', 'xb'):
            f = SNEWPYFlux(sn_model, flv, starttime, endtime)
            self.components[flv] = [f]

        return self


def get_starttime(starttime, minimum):
    """
    Return sensible start time.

    starttime - user-defined (defaults to None)
    minimum - earliest time in the input file
    """
    if starttime is None:
        starttime = ceil(minimum)
    elif starttime < minimum:
        raise ValueError("Start time cannot be earlier than %f (first entry in input file)." % minimum)
    return starttime


def get_endtime(endtime, maximum):
    """
    Return sensible end time.

    endtime - user-defined (defaults to None)
    maximum - latest time in the input file
    """
    if endtime is None:
        endtime = floor(maximum)
    elif endtime > maximum:
        raise ValueError("End time cannot be later than %f (last entry in input file)." % maximum)
    return endtime
