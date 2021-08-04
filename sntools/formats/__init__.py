from abc import ABC, abstractmethod
from importlib import import_module
from math import ceil, floor, pi


class BaseFlux(ABC):
    """Abstract base class that defines the interface for all submodules in `sntools.formats`."""
    def __repr__(self) -> str:
        if hasattr(self, '_repr'):
            return self._repr
        else:
            return super().__repr__()

    @abstractmethod
    def parse_input(self, input, inflv, starttime, endtime):
        """Read simulations data from input file.

        Arguments:
        input -- prefix of file containing neutrino fluxes
        inflv -- neutrino flavor to consider
        starttime -- start time set by user via command line option (or None)
        endtime -- end time set by user via command line option (or None)
        """
        # self.raw_times, self.starttime and self.endtime must be set in parse_input().
        # See get_raw_times(), get_starttime() and get_endtime() below.
        pass

    @abstractmethod
    def prepare_evt_gen(self, binned_t):
        """Pre-compute values necessary for event generation.

        Scipy/numpy are optimized for parallel operation on large arrays, making
        it orders of magnitude faster to pre-compute all values at one time
        instead of computing them lazily when needed.

        Argument:
        binned_t -- list of time bins for generating events
        """
        pass

    @abstractmethod
    def nu_emission(self, eNu, time):
        """Number of neutrinos emitted, as a function of energy.

        This does not include the geometry factor 1/(4 pi r**2).

        Arguments:
        eNu -- neutrino energy
        time -- time ;)
        """
        pass


class WeightedFlux(BaseFlux):
    def __init__(self, unweighted_flux, scale, distance=10.0) -> None:
        """Initialize a WeightedFlux.

        Arguments:
        unweighted_flux -- BaseFlux instance corresponding to original flux (at supernova)
        scale -- scaling factor from flux transformation, e.g. 1.0 for NoTransformation
        distance -- distance in kpc
        """
        super().__init__()
        self._repr = f"WeightedFlux({unweighted_flux}, scale={scale}, distance={distance})"
        self._flux = unweighted_flux
        distance *= 1.563738e32  # 1 kpc/(hbar * c) in MeV**(-1)
        self._weight = scale / (4 * pi * distance ** 2)

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
        self._repr = f"CompositeFlux.from_file('{file}', format='{format}', starttime={starttime}, endtime={endtime})"
        format = import_module('sntools.formats.' + format)

        for flv in ('e', 'eb', 'x', 'xb'):
            f = format.Flux()
            f._repr = repr(f)[:-1] + f": file='{file}', flavor='{flv}', starttime={starttime}, endtime={endtime}>"
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
                    wf = WeightedFlux(flux, scale, distance)
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
        self.raw_times = get_raw_times(times, self.starttime, self.endtime)

    def parse_input(self, *args):
        # handled in SNEWPYCompositeFlux.from_file()
        pass

    def prepare_evt_gen(self, *args):
        pass

    def nu_emission(self, eNu, time):
        from astropy import units as u
        nl = self._sn_model.get_initialspectra(time * u.ms, eNu * u.MeV, flavors=[self._flv])[self._flv]
        # SNEWPY uses cgs units internally, so convert before returning
        return nl.to(1 / u.MeV / u.ms).value


class SNEWPYCompositeFlux(CompositeFlux):
    """Adapter class to turn a SNEWPY.models.SupernovaModel into an sntools.formats.CompositeFlux"""

    @classmethod
    def from_file(cls, file, format, starttime=None, endtime=None):
        """Create a SNEWPYCompositeFlux from an input file."""
        self = SNEWPYCompositeFlux()
        self._repr = f"SNEWPYCompositeFlux.from_file('{file}', format='{format}', starttime={starttime}, endtime={endtime})"

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
        raise ValueError(f"Start time cannot be earlier than {minimum} (first entry in input file).")
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
        raise ValueError(f"End time cannot be later than {maximum} (last entry in input file).")
    return endtime


def get_raw_times(times, starttime, endtime):
    """
    Return list of times covering only the relevant time range.

    times - all times in input file
    starttime, endtime - relevant time range
    """
    i_min, i_max = 0, len(times) - 1
    for (i, time) in enumerate(times):
        if time < starttime:
            i_min = i
        elif time > endtime:
            i_max = i
            break
    return times[i_min: i_max + 1]
