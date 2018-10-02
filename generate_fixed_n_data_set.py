"""
Generate a fixed number of events (independent of supernova distance)
instead of a random (poisson-distributed) number of events for a fixed
SN distance (which is what `genevts.py` does).

This is designed to be a drop-in replacement for `genevts.py`; the only
difference is that `--distance` does not set the distance in kpc but
the number of events to generate.
"""
import __builtin__
from datetime import datetime
from importlib import import_module
import random
from scipy import integrate, interpolate

from channel import setup, ddEventRate, get_eNu, get_direction
from genevts import mixings, channels, detectors, parse_command_line_options


def main():
    args = parse_command_line_options()

    hierarchy = args.hierarchy
    global channels
    channels = channels if args.channel == "all" else [args.channel]
    input = args.input_file
    format = args.format
    output = args.output
    detector = detectors[args.detector]
    n_total = int(args.distance) # number of events to generate (see comment above)
    starttime = args.starttime if args.starttime else None
    endtime = args.endtime if args.endtime else None
    verbose = args.verbose


    channel_modules = {ic: import_module("interaction_channels." + ic) for ic in channels}

    # Take into account hierarchy-dependent flavor mixing and let channel.py
    # generate the actual events for each channel.
    binned_nevt_by_channel = {}
    for channel in channels:
        mod_channel = channel_modules[channel]
        for (original_flv, scale, detected_flv) in mixings[hierarchy]:
            if detected_flv in mod_channel.possible_flavors:

                if channel == "es": __builtin__._flavor = detected_flv

                # ignore dependence of `scale` on distance and detector size
                # since normalization is done by n_total instead

                cmd = "expected_nevts(_channel='%s', input='%s', _format='%s', inflv='%s', scale=%s, starttime=%s, endtime=%s, verbose=%s)" \
                    % (channel, input, format, original_flv, scale, starttime, endtime, verbose)
                if verbose: print "Now executing:", cmd
                binned_nevt_by_channel[(channel, original_flv, detected_flv)] = eval(cmd)


    events_to_generate = {}
    n_max = max([max(binned_events) for (binned_events, _, _, _) in binned_nevt_by_channel.values()])
    for _ in range(n_total):
        # pick a random channel & time bin, then accept it with a probability
        # proportional to the number of events in that time bin
        while True:
            key = random.choice(binned_nevt_by_channel.keys())
            (binned_nevt_th, _, _, _) = binned_nevt_by_channel[key]
            bin = random.randrange(len(binned_nevt_th))
            if n_max * random.random() < binned_nevt_th[bin]:
                break

        events_to_generate[key] = events_to_generate.get(key, []) + [bin]


    # iterate over `events_to_generate`
    events = []
    for (key, time_bins) in events_to_generate.iteritems():
        # TODO: setup stuff based on `key`
        (channel, original_flv, detected_flv) = key
        mod_channel = channel_modules[channel]

        # ensure that we're using fluxes for the right flavour
        _, binned_t, _starttime, _ = expected_nevts(channel, input, format, original_flv, scale, starttime, endtime, verbose)

        # generate the events
        for bin in time_bins:
            t = _starttime + (bin + random.random()) * bin_width
            eNu = get_eNu(binned_t[bin])
            (dirx, diry, dirz) = get_direction(eNu)
            eE = mod_channel.get_eE(eNu, dirz)
            events.append((t, mod_channel.pid, eE, dirx, diry, dirz))

    with open(output, 'w') as outfile:
        outfile.write("# Generated on %s with the options:\n" % datetime.now())
        outfile.write("# " + str(args) + "\n")
        for (t, pid, e, dirx, diry, dirz) in events:
            outfile.write("%.5f, %i, %.5f, %.5f, %.5f, %.5f\n" % (t, pid, e, dirx, diry, dirz))


def expected_nevts(_channel, input, _format, inflv, scale, starttime, endtime, verbose):
    """Calculate expected number of events by interpolating from time steps in the input data.

    Arguments:
    _channel -- abbreviation of interaction channel, e.g. 'ibd', 'es', ...
    input -- name (or common prefix) of file(s) containing neutrino fluxes
    _format -- which parser (in folder `formats/`) to use for input file(s)
    inflv -- original neutrino flavor (at time of production in the SN)
    scale -- constant factor, accounts for oscillation probability, distance of SN, size of detector
    starttime -- start time set by user via command line option (or None)
    endtime -- end time set by user via command line option (or None)
    """
    channel, format = setup(_channel, _format) # import appropriate modules
    scale *= channel.targets_per_molecule

    (starttime, endtime, raw_times) = format.parse_input(input, inflv, starttime, endtime)

    # integrate over eE and then eNu to obtain the event rate at time t
    raw_nevts = [scale * integrate.nquad(ddEventRate, [channel.bounds_eE, channel.bounds_eNu], args=[t])[0]
                 for t in raw_times]
    event_rate = interpolate.pchip(raw_times, raw_nevts)

    global bin_width
    bin_width = 1 # in ms
    n_bins = int((endtime - starttime)/bin_width) # number of full-width bins; int() implies floor()

    # scipy is optimized for operating on large arrays, making it orders of
    # magnitude faster to pre-compute all values of the interpolated functions.
    binned_t = [starttime + (i+0.5)*bin_width for i in range(n_bins)]
    binned_nevt_th = event_rate(binned_t)
    format.prepare_evt_gen(binned_t) # give flux script a chance to pre-compute values

    return binned_nevt_th, binned_t, starttime, endtime


if __name__ == "__main__":
    main()
