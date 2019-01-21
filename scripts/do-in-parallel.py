import argparse
from multiprocessing import Pool
import os
import subprocess
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("step", help="What step of the analysis chain to do. \
                                  Can be 'generate_dataset', 'reconstruct_energy' or 'readd_time'.")
args = parser.parse_args()

# input files: one per model family (to study explosion mechanism)
inputs = (
          ("in_tamborra27", "garching"),
          ("in_spectob2000.data", "nakazato"),
          ("in_totani20", "totani"),
         )

nevts = (3e5,) #(300, 1000, 3000, 10000, 30000)

hierarchies = ("noosc", "normal", "inverted")

imax = 1

def generate_dataset(input, nevt, hierarchy, i):
    print "... now running:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.split('.')[0] # remove file extension, if any
    path = "out/%s/%s" % (nevt, format)
    outfile = "%s/%s-%s-%s.csv" % (path, _infile, hierarchy, i)

    if not os.path.isdir(path):
        os.makedirs(path)
    elif os.path.isfile(outfile):
        print('File %s already exists. Aborting to avoid overwriting it ...' % outfile)
        raise Exception


    cmd = "python generate_fixed_n_data_set.py %s -f %s -o %s --hierarchy=%s --distance %i \
           --starttime=20 --endtime=520" % (infile, format, outfile, hierarchy, nevt)
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


def reconstruct_energy(input, nevt, hierarchy, i):
    print "... now reconstructing:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.split('.')[0] # remove file extension, if any
    path = "out/%s/%s" % (nevt, format)
    outfile = "%s/%s-%s-%s" % (path, _infile, hierarchy, i)

    if not os.path.isdir(path):
        print("Directory %s doesn't exist. Aborting ..." % path)
        raise Exception
    elif os.path.isfile(outfile + ".reco"):
        print('File %s.reco already exists. Aborting to avoid overwriting it ...' % outfile)
        raise Exception

    cmd = "root -q '../energetic-bonsai/energetic_bonsai_final.C(\"%s.root\")' > %s.reco" % (outfile, outfile)
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


def readd_time(input, nevt, hierarchy, i):
    infile, format = input
    _infile = infile.split('.')[0] # remove file extension, if any
    path = "out/%s/%s" % (nevt, format)
    outfile = "%s/%s-%s-%s" % (path, _infile, hierarchy, i)

    if not os.path.isfile(outfile + '.reco'):
        print("File %s.reco doesn't exist. Aborting ..." % outfile)
        raise Exception
    elif os.path.isfile(outfile + '.reco2'):
        print('File %s.reco2 already exists. Aborting to avoid overwriting it ...' % outfile)
        raise Exception


    cmd = "python ../sntools/scripts/add-time-to-reconstructed-evts.py %s" % outfile
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


if args.step == 'generate_dataset':
    f = generate_dataset
elif args.step == 'reconstruct_energy':
    f = reconstruct_energy
elif args.step == 'readd_time':
    f = readd_time
else:
    raise Exception("%s: No such function." % args.step)


pool = Pool(processes=4) # start 4 worker processes

results = []
for input in inputs:
    for nevt in nevts:
        for hierarchy in hierarchies:
            for i in range(1, imax+1):
                print "Adding to pool:", input, str(nevt), hierarchy, str(i)
                res = pool.apply_async(f, (input, nevt, hierarchy, i))
                results.append(res)
                sleep(1) # this appears to fix an issue where the pool would skip some of the first tasks

for res in results:
#     print res.get(timeout=60)
    res.wait()
