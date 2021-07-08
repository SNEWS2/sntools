#!/usr/bin/python

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from importlib import import_module
import random

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


def main():
    args = parse_command_line_options()
    random.seed(args.randomseed)

    # Calculate fluxes at detector
    raw_flux = CompositeFlux.from_file(args.input_file, args.format, args.starttime, args.endtime)
    flux_at_detector = raw_flux.transformed_by(args.transformation, args.distance)

    # Generate events for each (sub-)channel and combine them
    pool = ProcessPoolExecutor(max_workers=args.maxworkers)
    results = []
    for channel in sorted(args.channels):
        mod_channel = import_module("sntools.interaction_channels." + channel)
        n_targets = args.detector.n_molecules * args.detector.material["channel_weights"][channel]
        for flv in mod_channel.possible_flavors:
            channel_instance = mod_channel.Channel(flv)
            for flux in flux_at_detector.components[flv]:
                results.append(pool.submit(gen_evts, channel_instance, flux, n_targets, args.randomseed + random.random(), args.verbose))

    events = []
    for result in as_completed(results):
        events.extend(result.result())

    # Sort events by time and write them to an output file
    events.sort(key=lambda evt: evt.time)
    with open(args.output, "w") as outfile:
        if args.verbose:  # write parameters to file as a comment
            outfile.write("# Generated on %s with the options:\n" % datetime.now())
            outfile.write("# " + str(args) + "\n")
        if args.mcformat == 'NUANCE':
            for (i, evt) in enumerate(events):
                evt.vertex = args.detector.generate_random_vertex()
                outfile.write(evt.nuance_string(i))
            outfile.write("$ stop\n")
        if args.mcformat == 'RATPAC':
            for (i, evt) in enumerate(events):
                evt.vertex = args.detector.generate_random_vertex()
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

    # NOTE: Deprecated since v1.0 in favor of '--transformation'
    choices = ("noosc", "normal", "inverted")
    parser.add_argument("-H", "--hierarchy", "--ordering", choices=choices, help=argparse.SUPPRESS, action=DeprecationAction)

    choices = ("NoTransformation", "AdiabaticMSW_NMO", "AdiabaticMSW_IMO")
    parser.add_argument("-t", "--transformation", metavar="TRANSFORM", choices=choices, default=choices[0],
                        help="Transformation between neutrino flux inside SN and flux in the detector on Earth. \
                              Choices: %(choices)s. Default: %(default)s.")

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

    args = parser.parse_args()

    args.transformation = Transformation(args.transformation)
    args.detector = Detector(args.detector)
    args.channels = args.detector.material["channel_weights"] if args.channel == "all" else [args.channel]
    del args.hierarchy  # see args.transformation
    del args.channel  # see args.channels

    if args.verbose:
        print("Arguments:")
        for (k, v) in vars(args).items():
            print(f"  {k}: {v}")

    return args


class DeprecationAction(argparse.Action):
    """argparse.Action subclass to deprecate options"""

    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if self.dest == 'hierarchy':
            replacement = {"noosc": "NoTransformation", "normal": "AdiabaticMSW_NMO", "inverted": "AdiabaticMSW_IMO"}[values]
            print(f"❌ '{option_string} {values}' is deprecated. Please switch to using '--transformation {replacement}'.")
            self.dest = 'transformation'
            values = replacement
        else:
            print(f"❌ '{option_string}' is deprecated. No replacement available.")

        setattr(namespace, self.dest, values)


if __name__ == "__main__":
    main()
