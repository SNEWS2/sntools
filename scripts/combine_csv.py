"""
Combines multiple csv files (one for each interaction channel / input flavour
combination), each containing event counts for multiple time and energy bins,
into one file.
This allows creating one file for each hierarchy (noosc/normal/inverted) and
then perform a chi-squared analysis on that instead of on a single channel.
"""
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Name of input directory containing 1 csv file per interaction channel.")
parser.add_argument("output", help="Name of output csv file.")
args = parser.parse_args()


inputs = [os.path.join(args.input, file) for file in os.listdir(args.input)
                                         if os.path.splitext(file)[1] in (".csv", ".txt")]

ntot_dict, nevt_dict = {}, {}
first_time = True

for input in inputs:
    # read from input file
    timebins = []
    with open(input) as infile:
        for line in infile:
            if (line.startswith("#") or line.isspace()):
                comment = line # note on format is the last comment line
            else:
                timebins.append(map(float, line.split(",")))

    for timebin in timebins:
        t, n_tot, n_binned = timebin[0], timebin[1], timebin[2:]
        if first_time:
            ntot_dict[t] = n_tot
            nevt_dict[t] = n_binned
        else:
            ntot_dict[t] += n_tot
            for i, new in enumerate(n_binned):
                nevt_dict[t][i] += new

    first_time = False

with open(args.output, 'w') as outfile:
    outfile.write(comment + '\n')
    for t in sorted(nevt_dict):
        outfile.write(str(t) + ', ' + str(ntot_dict[t]))
        for n in nevt_dict[t]:
            outfile.write(', ' + str(n))
        outfile.write('\n')
