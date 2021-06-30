from abc import ABC, abstractmethod
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
