#!/usr/bin/python

from __future__ import print_function

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from importlib import import_module
import random

import numpy as np

try:
    import sntools  # if sntools was installed via pip
except ImportError:
    # if running this directly from the repo, modify `sys.path` to ensure all imports work
    import os
    import sys
    abs_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(abs_dir)
    sys.path.append(parent_dir)
    import sntools

from sntools.channel import gen_evts
from sntools.detectors import Detector, supported_detectors
from sntools.formats import CompositeFlux
from sntools.transformation import Transformation

import sys
if sys.version_info < (3, 6):
    from sntools import tryprint
    tryprint(u"\u274c", "[WARNING]")
    print("You are using Python %s.%s.%s, which is not supported any more." %
          (sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    tryprint(u"\u274c", "[WARNING]")
    print("Future versions of sntools will require Python 3.6 or higher and it is recommended that you transition as soon as possible.")
    tryprint(u"\u274c", "[WARNING]")
    print("See https://github.com/JostMigenda/sntools/issues/29 for more information.\n")


def main():
    args = parse_command_line_options()

    detector = Detector(args.detector)
    channels = detector.material["channel_weights"] if args.channel == "all" else [args.channel]
    transformation = Transformation(args.hierarchy)
    input = args.input_file
    format = args.format
    output = args.output
    mcformat = args.mcformat
    distance = args.distance
    starttime = args.starttime
    endtime = args.endtime
    seed = args.randomseed
    verbose = args.verbose

    if verbose:
        print("channel(s) =", channels)
        print("transform. =", transformation)
        print("input file =", input, "--- format =", format)
        print("output     =", output)
        print("mcformat   =", mcformat)
        print("detector   =", detector)
        print("distance   =", distance)
        print("starttime  =", starttime)
        print("endtime    =", endtime)
        print("randomseed =", seed)
        print("**************************************")

    random.seed(seed)

    # Calculate fluxes at detector
    raw_flux = CompositeFlux.from_file(args.input_file, args.format, args.starttime, args.endtime)
    flux_at_detector = raw_flux.transformed_by(transformation, args.distance)

    # Generate events for each (sub-)channel and combine them
    pool = ProcessPoolExecutor(max_workers=args.maxworkers)
    results = []
    for channel in sorted(channels):
        mod_channel = import_module("sntools.interaction_channels." + channel)
        n_targets = detector.n_molecules * detector.material["channel_weights"][channel]
        for flv in mod_channel.possible_flavors:
            channel_instance = mod_channel.Channel(flv)
            for flux in flux_at_detector.components[flv]:
                results.append(pool.submit(gen_evts, channel_instance, flux, n_targets, seed + random.random(), args.verbose))

    events = []
    for result in as_completed(results):
        events.extend(result.result())

    # Sort events by time and write them to a nuance-formatted output file
    events.sort(key=lambda evt: evt.time)
    with open(output, "w") as outfile:
        if verbose:  # write parameters to file as a comment
            outfile.write("# Generated on %s with the options:\n" % datetime.now())
            outfile.write("# " + str(args) + "\n")
        if mcformat == 'NUANCE':
            for (i, evt) in enumerate(events):
                evt.vertex = detector.generate_random_vertex()
                outfile.write(evt.nuance_string(i))
            outfile.write("$ stop\n")
        if mcformat == 'RATPAC':
            for (i, evt) in enumerate(events):
                evt.vertex = detector.generate_random_vertex()
                outfile.write(evt.ratpac_string(i, events))


def parse_command_line_options():
    """Define and parse command line options."""
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help="Name or common prefix of the input file(s). Required.")

    choices = ["gamma", "nakazato", "princeton", "totani", "warren2020"]
    default = "totani"
    parser.add_argument("-f", "--format", metavar="FORMAT", choices=choices, default=default,
                        help="Format of input files. See parsers in folder 'formats/' \
                              for details. Choices: %s. Default: '%s'." % (choices, default))

    default = "outfile.kin"
    parser.add_argument("-o", "--output", metavar="FILE", default=default,
                        help="Name of the output file. Default: '%s'." % default)

    choices = ["NUANCE", "RATPAC"]
    default = "NUANCE"
    parser.add_argument("-m", "--mcformat", metavar="MCFORMAT", choices=choices, default=default,
                        help="MC output format for simulations. Choices: %s. Default: %s." % (choices, default))

    choices = ("noosc", "normal", "inverted")
    parser.add_argument("-H", "--hierarchy", "--ordering", choices=choices, help=argparse.SUPPRESS)

    choices = ["ibd", "es", "o16e", "o16eb", "c12e", "c12eb", "c12nc"]
    parser.add_argument("-c", "--channel", metavar="INTCHANNEL", choices=choices, default="all",
                        help="Interaction channels to consider. Currently, inverse beta decay (ibd), \
                              electron scattering (es), nu_e + oxygen CC (o16e), nu_e-bar + oxygen CC \
                              (o16eb), nu_e + carbon CC (c12e), nu_e-bar + carbon CC (c12eb) \
                              and nu + carbon NC (c12nc) are supported. \
                              Choices: %s. Default: All supported channels." % choices)

    default = "HyperK"
    parser.add_argument("-d", "--detector", metavar="DETECTOR", choices=supported_detectors, default=default,
                        help="Detector configuration. Choices: %s. Default: '%s'." % (supported_detectors, default))

    default = 10.0
    parser.add_argument("--distance", type=float, default=default,
                        help="Distance to supernova in kpc. Default: '%s'." % default)

    parser.add_argument("--starttime", metavar="T", type=float,
                        help="Start generating events at T milliseconds. Default: First time bin in input file.")

    parser.add_argument("--endtime", metavar="T", type=float,
                        help="Stop generating events at T milliseconds. Default: Last time bin in input file.")

    parser.add_argument("--randomseed", metavar="SEED", type=int,  # non-ints may not give reproducible results
                        help="Integer used as a random number seed to reproducibly generate events. Default: None.")

    parser.add_argument("--maxworkers", metavar="N", type=int,
                        help="Maximum number of parallel processes. Default: [number of CPU cores].")

    parser.add_argument("-v", "--verbose", action="count", help="Verbose output, e.g. for debugging. Off by default.")

    parser.add_argument('--version', action='version', version='%(prog)s ' + sntools.__version__)

    return parser.parse_args()


if __name__ == "__main__":
    main()
