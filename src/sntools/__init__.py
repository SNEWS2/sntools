"""An Event Generator for Supernova Burst Neutrinos

sntools is a Monte Carlo event generator for supernova and pre-supernova neutrino interactions.
Based on detailed time- and energy-dependent neutrino fluxes provided by
various supernova models, it generates interactions within the detector volume
and writes them to event files that can be used as an input for a full detector
simulation.
sntools was originally developed for Hyper-Kamiokande using core collapse supernova models and later extended to 
support different detectors, detector materials and to use pre supernova models.

For usage information, run `sntools -h` or `python sntools/genevts.py -h`.
For more extensive documentation, to report issues or to contribute code,
see https://github.com/SNEWS2/sntools.
"""
__version__ = '1.5b1'


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
    flux_url = 'https://raw.githubusercontent.com/SNEWS2/sntools/main/fluxes/intp2001.data'
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

    test_sha = "0e3444d5a174ee49bd3adad498a2a022d2340736f2e056bfab12df0e5a9b0c54"
    if output_sha == test_sha:
        tryprint(u"\u2705", "[SUCCESS]")
        print("Everything seems to work fine. Enjoy using sntools!")
    else:
        tryprint(u"\u274c", "[ERROR]")
        print("Test did not generate the expected events.")
        tryprint(u"\u274c", "[ERROR]")
        print("Please ensure you have installed the most recent version of sntools and all dependencies.")
        tryprint(u"\u274c", "[ERROR]")
        print("If this persists, please go to https://github.com/SNEWS2/sntools and open a new issue.")

def presnsetup():
    """
    Downloads sample flux file from GitHub if necessary and performs pre supernova integration test.
    """
    tryprint(u"\u2705", "[SUCCESS]")
    print("sntools was imported from " + __path__[0])
    import hashlib
    import os
    import sys
    from . import genevts

    flux_dir = 'fluxes/'
    flux_file = flux_dir + 'totalLuminosity_15SolarMass.dat'
    flux_url = 'https://raw.githubusercontent.com/SNEWS2/snewpy-models-presn/master/models/Patton_2017/totalLuminosity_15SolarMass.dat'
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
    sys.argv += [flux_file, '--format', 'SNEWPY-Patton_2017', '--detector', 'HyperK', '--distance', '0.15', '--transformation', 'AdiabaticMSW_IMO', '--starttime', '-120', '--endtime', '0', '--binsize', '1', '-o', 'presnoutfile.kin', '--randomseed', '205']
    genevts.main()

    tryprint(u"\U0001f6e0")
    print("Checking output file ...")
    with open('presnoutfile.kin', 'r') as f:
        output_sha = hashlib.sha256(f.read().encode('utf-8')).hexdigest()

    test_sha = "824dffa7ed01016a9044a28b36ca33aacdd81a9db662ca69a08d9ab0eb32bf29"
    if output_sha == test_sha:
        tryprint(u"\u2705", "[SUCCESS]")
        print("Everything seems to work fine. Enjoy using sntools!")
    else:
        tryprint(u"\u274c", "[ERROR]")
        print("Test did not generate the expected events.")
        tryprint(u"\u274c", "[ERROR]")
        print("Please ensure you have installed the most recent version of sntools and all dependencies.")
        tryprint(u"\u274c", "[ERROR]")
        print("If this persists, please go to https://github.com/SNEWS2/sntools and open a new issue.")

def presnbinsizetest():
    """
    Downloads sample flux file from GitHub if necessary and performs pre supernova binsize test.
    """
    tryprint(u"\u2705", "[SUCCESS]")
    print("sntools was imported from " + __path__[0])
    import hashlib
    import os
    import sys
    from . import genevts

    flux_dir = 'fluxes/'
    flux_file = flux_dir + 'totalLuminosity_15SolarMass.dat'
    flux_url = 'https://raw.githubusercontent.com/SNEWS2/snewpy-models-presn/master/models/Patton_2017/totalLuminosity_15SolarMass.dat'
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
    sys.argv += [flux_file, '--format', 'SNEWPY-Patton_2017', '--detector', 'HyperK', '--distance', '0.15', '--transformation', 'AdiabaticMSW_NMO', '--starttime', '-120', '--endtime', '0', '--binsize', '10', '-o', 'presnbinsizeoutfile.kin', '--randomseed', '100']
    genevts.main()

    tryprint(u"\U0001f6e0")
    print("Checking output file ...")
    with open('presnbinsizeoutfile.kin', 'r') as f:
        output_sha = hashlib.sha256(f.read().encode('utf-8')).hexdigest()

    test_sha = "7e1e8d14c732d171a9ad9f924bde4ca31f6897ca24aa0c94074e008fd69e76e7"
    if output_sha == test_sha:
        tryprint(u"\u2705", "[SUCCESS]")
        print("Everything seems to work fine. Enjoy using sntools!")
    else:
        tryprint(u"\u274c", "[ERROR]")
        print("Test did not generate the expected events.")
        tryprint(u"\u274c", "[ERROR]")
        print("Please ensure you have installed the most recent version of sntools and all dependencies.")
        tryprint(u"\u274c", "[ERROR]")
        print("If this persists, please go to https://github.com/SNEWS2/sntools and open a new issue.")


def tryprint(default, alternative=''):
    try:
        print(default, end=' ')
    except UnicodeEncodeError:
        print(alternative, end=' ')
