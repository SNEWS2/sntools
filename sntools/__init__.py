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

__version__ = '0.7.2.post2'


def setup():
    """
    Downloads sample flux file from GitHub if necessary and performs integration test.
    """
    print(u"\u2705 sntools was imported from " + __path__[0])
    import hashlib
    import os
    import sys
    from . import genevts

    flux_dir = 'fluxes/'
    flux_file = flux_dir + 'intp2001.data'
    flux_url = 'https://raw.githubusercontent.com/JostMigenda/sntools/master/fluxes/intp2001.data'
    if os.path.exists(flux_file):
        print(u"\u2705 Using sample flux file at " + flux_file)
    else:
        print(u"\U0001f6e0 Downloading sample flux file from " + flux_url)
        if not os.path.isdir(flux_dir):
            os.mkdir(flux_dir)
        try:
            from urllib.request import urlretrieve  # Python 3.x
        except ImportError:
            from urllib import urlretrieve  # Python 2.7
        try:
            urlretrieve(flux_url, filename=flux_file)
            print(u"\u2705 Saved sample flux file to " + flux_file)
        except IOError:
            print(u"\u274c Error: Cannot download sample flux file.")
            sys.exit(-1)

    print(u"\U0001f6e0 Testing event generation ...")
    sys.argv += [flux_file, '--format', 'nakazato', '--detector', 'WATCHMAN-LS', '--distance', '2', '--ordering', 'normal', '--starttime', '100', '--endtime', '300', '-o', 'outfile.kin', '--randomseed', '314']
    genevts.main()

    print(u"\U0001f6e0 Checking output file ...")
    with open('outfile.kin', 'rb') as f:
        output_sha = hashlib.sha256(f.read()).hexdigest()
    
    test_sha = "6173ea086aab0590212d0207612d57a46184f7569d0205d322fe3b5689f5f700"
    if output_sha == test_sha:
        print(u"\u2705 Everything seems to work fine. Enjoy using sntools!")
    else:
        print(u"\u274c Error: Test did not generate the expected events.\n"
              u"\u274c Please ensure you have installed the most recent version of sntools and all dependencies.\n"
              u"\u274c If this persists, please go to https://github.com/JostMigenda/sntools and open a new issue.")
