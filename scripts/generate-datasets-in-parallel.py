from multiprocessing import Pool
import os
import subprocess
from time import sleep

# TODO: add more files: Garching directions, Nakazato masses/metallicities
inputs = (("infile_full", "garching"),
#           ("infile_intp2002.txt", "nakazato"),
          ("infile_wilson", "totani")
         )

nevts = (300, 1000, 3000)#, 10000, 30000)

hierarchies = ("noosc", "normal", "inverted")

imax = 5

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

pool = Pool(processes=4) # start 4 worker processes

results = []
for input in inputs:
    for nevt in nevts:
        for hierarchy in hierarchies:
            for i in range(1, imax+1):
                print "Adding to pool:", input, str(nevt), hierarchy, str(i)
                res = pool.apply_async(generate_dataset, (input, nevt, hierarchy, i))
                results.append(res)
                sleep(1) # this appears to fix an issue where the pool would skip some of the first tasks

for res in results:
#     print res.get(timeout=60)
    res.wait()
