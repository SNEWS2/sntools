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
__version__ = '1.0b2'


def setup():
    """
    Downloads sample flux file from GitHub if necessary and performs integration test.
    """
    tryprint(u"\u2705", "[SUCCESS]")
    print("sntools was imported from " + __path__[0])
    import hashlib
    import os
    import sys
    from . import genevts

    flux_dir = 'fluxes/'
    flux_file = flux_dir + 'intp2001.data'
    flux_url = 'https://raw.githubusercontent.com/JostMigenda/sntools/master/fluxes/intp2001.data'
    if os.path.exists(flux_file):
        tryprint(u"\u2705", "[SUCCESS]")
        print("Using sample flux file at " + flux_file)
    else:
        tryprint(u"\U0001f6e0")
        print("Downloading sample flux file from " + flux_url)
        if not os.path.isdir(flux_dir):
            os.mkdir(flux_dir)

        from urllib.request import urlretrieve
        try:
            urlretrieve(flux_url, filename=flux_file)
            tryprint(u"\u2705", "[SUCCESS]")
            print("Saved sample flux file to " + flux_file)
        except IOError:
            tryprint(u"\u274c", "[ERROR]")
            print("Cannot download sample flux file.")
            sys.exit(-1)

    tryprint(u"\U0001f6e0")
    print("Testing event generation ...")
    sys.argv += [flux_file, '--format', 'nakazato', '--detector', 'WATCHMAN-LS', '--distance', '2', '--transformation', 'AdiabaticMSW_NMO', '--starttime', '100', '--endtime', '300', '-o', 'outfile.kin', '--randomseed', '314']
    genevts.main()

    tryprint(u"\U0001f6e0")
    print("Checking output file ...")
    with open('outfile.kin', 'r') as f:
        output_sha = hashlib.sha256(f.read().encode('utf-8')).hexdigest()

    test_sha = "4424b738d32b6c4baa459fb6fb76d63f944b572a02d319bde4fd247a520f1e29"
    if output_sha == test_sha:
        tryprint(u"\u2705", "[SUCCESS]")
        print("Everything seems to work fine. Enjoy using sntools!")
    else:
        tryprint(u"\u274c", "[ERROR]")
        print("Test did not generate the expected events.")
        tryprint(u"\u274c", "[ERROR]")
        print("Please ensure you have installed the most recent version of sntools and all dependencies.")
        tryprint(u"\u274c", "[ERROR]")
        print("If this persists, please go to https://github.com/JostMigenda/sntools and open a new issue.")


def tryprint(default, alternative=''):
    try:
        print(default, end=' ')
    except UnicodeEncodeError:
        print(alternative, end=' ')
