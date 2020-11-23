#!/usr/bin/python

from __future__ import print_function

try:
    import __builtin__ as builtins  # Python 2.7
except ImportError:
    import builtins  # Python 3
import argparse
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


# mixing parameters from P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020, 083C01 (2020)
s12 = 0.307  # sin^2 theta_12
c12 = 1 - s12
s13 = 0.0218  # sin^2 theta_13
c13 = 1 - s13

# While exiting the supernova, neutrinos experience mass hierarchy-dependent
# flavor transitions via the MSW effect. This dictionary contains 3-tuples of
#   * original flavor at production (i.e. in input files from computer simulations),
#   * mixing probability,
#   * resulting flavor in detector.
# See the section "Treatment of Neutrino Flavour Conversion" in the documentation.
mixings = {
    "noosc": (
        ("e", 1, "e"),
        ("eb", 1, "eb"),
        ("x", 2, "x"),  # scale = 2 to include both nu_mu and nu_tau
        ("xb", 2, "xb"),
    ),
    "normal": (
        ("e", s13, "e"),  # nu_e that originated as nu_e
        ("x", c13, "e"),  # nu_e that originated as nu_x
        ("eb", c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_e
        ("xb", 1 - c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_x
        ("e", c13, "x"),  # nu_x that originated as nu_e
        ("x", 1 + s13, "x"),  # nu_x that originated as nu_x
        ("eb", 1 - c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_e
        ("xb", 1 + c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_x
    ),
    "inverted": (
        ("e", s12 * c13, "e"),  # nu_e that originated as nu_e
        ("x", 1 - s12 * c13, "e"),  # nu_e that originated as nu_x
        ("eb", s13, "eb"),  # anti-nu_e that originated as anti-nu_e
        ("xb", c13, "eb"),  # anti-nu_e that originated as anti-nu_x
        ("e", 1 - s12 * c13, "x"),  # nu_x that originated as nu_e
        ("x", 1 + s12 * c13, "x"),  # nu_x that originated as nu_x
        ("eb", c13, "xb"),  # anti-nu_x that originated as anti-nu_e
        ("xb", 1 + s13, "xb"),  # anti-nu_x that originated as anti-nu_x
    ),
}


def main():
    args = parse_command_line_options()

    detector = Detector(args.detector)
    channels = detector.material["channel_weights"] if args.channel == "all" else [args.channel]
    hierarchy = args.hierarchy
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
        print("hierarchy  =", hierarchy)
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
    np.random.seed(seed)

    # Take into account hierarchy-dependent flavor mixing and let channel.py
    # generate the actual events for each channel.
    events = []
    for channel in sorted(channels):
        mod_channel = import_module("sntools.interaction_channels." + channel)
        for (original_flv, scale, detected_flv) in mixings[hierarchy]:
            if detected_flv in mod_channel.possible_flavors:

                # TODO: Replace this with a more sensible design, e.g. see https://stackoverflow.com/a/15959638
                builtins._flavor = detected_flv

                scale *= (10.0 / distance) ** 2  # flux is proportional to 1/distance**2
                scale *= detector.n_molecules
                scale *= detector.material["channel_weights"][channel]

                cmd = "gen_evts(_channel='%s', input='%s', _format='%s', inflv='%s', scale=%s, starttime=%s, endtime=%s, verbose=%s)" \
                    % (channel, input, format, original_flv, scale, starttime, endtime, verbose)
                if verbose:
                    print("Now executing:", cmd)
                events.extend(eval(cmd))

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

    choices = ["noosc", "normal", "inverted"]
    default = choices[0]
    parser.add_argument("-H", "--hierarchy", "--ordering", metavar="HIERARCHY", choices=choices, default=default,
                        help="Oscillation scenario. Choices: %s. Default: '%s'." % (choices, default))

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

    parser.add_argument("-v", "--verbose", action="count", help="Verbose output, e.g. for debugging. Off by default.")

    parser.add_argument('--version', action='version', version='%(prog)s ' + sntools.__version__)

    return parser.parse_args()


if __name__ == "__main__":
    main()
