import argparse
from multiprocessing import Pool
import os
import subprocess
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("step", help="What step of the analysis chain to do. Options: \
                                  1) 'generate_dataset' \
                                  2) 'split_eventfile' \
                                  3) 'csv_to_nuance' \
                                  4) 'simulate_events' \
                                  5) 'reconstruct_energy' \
                                  6) 'readd_time'.")
args = parser.parse_args()

# input files: one per model family (to study explosion mechanism)
inputs = (
          ("in_tamborra27.txt", "garching"),
          ("in_couch_a0.8_m20.0.dat", "garching"),
          ("in_spectob2000a.data", "nakazato"),
          #("in_spectob2000.data", "nakazato"),
          ("in_totani20", "totani"),
          ("in_vartanyan9.dat", "princeton"),
         )

# input files of different models in same model family (to study progenitor dependence)
# inputs = (
#           ("in_intp2002.data", "nakazato"),
#           ("in_intp2012.data", "nakazato"),
#           ("in_intp3002.data", "nakazato"),
#           ("in_intp3012.data", "nakazato"),
#          )

hierarchies = ("noosc", "normal", "inverted")

nevts = (100, 300)
imax = 1000
ntotal = max(nevts) * imax
detector = 'HyperK' # SuperK, HyperK

def generate_dataset(input, nevt, hierarchy, i):
    print "... now generating events:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)
    outfile = "%s/%s.csv" % (path, i)

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


def split_eventfile(input, nevt, hierarchy, i):
    print "... now splitting into multiple files:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any

    # input file with ntotal generated events
    inpath = "out/%s/%s/%s/1.csv" % (ntotal, _infile, hierarchy)
    # output files
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)

    if not os.path.isdir(path):
        os.makedirs(path)
    elif os.path.isfile(path + '/1.csv'):
        print('File %s/1.csv already exists. Aborting to avoid overwriting it ...' % path)
        raise Exception

    nfiles = ntotal / max(nevts)
    cmd = "python ../sntools/scripts/split_eventfile.py %s -o %s --nevts %s --nfiles %i" % (inpath, path, nevt, nfiles)
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


def csv_to_nuance(input, nevt, hierarchy, i):
    print "... now generating NUANCE files:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)

    if not os.path.isdir(path):
        print("Directory %s doesn't exist. Aborting ..." % path)
        raise Exception
    elif os.path.isfile(path + '/1.kin'):
        print('File %s/1.kin already exists. Aborting to avoid overwriting it ...' % path)
        raise Exception

    cmd = "python ../sntools/csv2nuance.py %s --wcsim %s/wcsim.mac --detector %s" % (path, path, detector)
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


def simulate_events(input, nevt, hierarchy, i):
    print "... now simulating:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)

    if not os.path.isdir(path):
        print("Directory %s doesn't exist. Aborting ..." % path)
        raise Exception
    elif os.path.isfile(path + "/1.root"):
        print('File %s/1.root already exists. Aborting to avoid overwriting it ...' % path)
        raise Exception

    cmd = "WCSim %s/wcsim.mac" % (path)
    print cmd
    result = subprocess.call(cmd, shell=True)
    return result


def reconstruct_energy(input, nevt, hierarchy, i):
    print "... now reconstructing:", input, str(nevt), hierarchy, str(i)
    infile, format = input
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)
    outfile = "%s/%s" % (path, i)

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
    _infile = infile.rsplit('.',1)[0] # remove file extension, if any
    path = "out/%s/%s/%s" % (nevt, _infile, hierarchy)
    outfile = "%s/%s" % (path, i)

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
    nevts = (ntotal,)
    imax = 1
elif args.step == 'split_eventfile':
    f = split_eventfile
    imax = 1
elif args.step == 'csv_to_nuance':
    f = csv_to_nuance
    imax = 1
elif args.step == 'simulate_events':
    f = simulate_events
    imax = 1
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
                sleep(0.5) # this appears to fix an issue where the pool would skip some of the first tasks

for res in results:
    res.wait()
