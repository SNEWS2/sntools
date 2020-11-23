"""An Event Generator for Supernova Burst Neutrinos

sntools is a Monte Carlo event generator for supernova neutrino interactions.
Based on detailed time- and energy-dependent neutrino fluxes provided by
various supernova models, it generates interactions within the detector volume
and writes them to event files that can be used as an input for a full detector
simulation.
sntools was originally developed for Hyper-Kamiokande and later extended to
support different detectors and detector materials.

For usage information, run `sntools -h` or `python sntools/genevts.py -h`.
For more extensive documentation, to report issues or to contribute code,
see https://github.com/JostMigenda/sntools.
"""

__version__ = '0.7.0'
