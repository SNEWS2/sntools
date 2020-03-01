import unittest

from . import gamma as f


class ParseInputTest(unittest.TestCase):
    def test_default_startendtime(self):
        # test extracting default start- and endtime from sample input file
        starttime, endtime, _ = f.parse_input('fluxes/sample-gamma.txt', 'e', None, None)
        self.assertEqual(starttime, 101)
        self.assertEqual(endtime, 125)

    def test_valid_startendtime(self):
        # test extracting default start- and endtime from sample input file
        starttime, endtime, times = f.parse_input('fluxes/sample-gamma.txt', 'e', 108, 115)
        self.assertEqual(starttime, 108.0)
        self.assertEqual(endtime, 115.0)
        self.assertEqual([round(t, 3) for t in times], # avoid floating point inaccuracies
                         [107.505, 108.001, 108.505, 109.002, 109.495, 110.005,
                          110.504, 111.003, 111.500, 112.395, 112.499, 113.004,
                          113.503, 114.002, 114.501, 114.999, 115.499])

    def test_invalid_startendtime(self):
        # should raise ValueError if specified start-/endtime is invalid
        with self.assertRaises(ValueError):
            f.parse_input('fluxes/sample-gamma.txt', 'e', 100, None)

        with self.assertRaises(ValueError):
            f.parse_input('fluxes/sample-gamma.txt', 'e', None, 130)

    def test_invalid_filename(self):
        # should raise IOError (FileNotFoundError in Python 3) if file doesn't exist
        with self.assertRaises(IOError):
            f.parse_input('fluxes/NONEXISTENT_FILE.data', 'e', None, None)


class PrepareEvtGenTest(unittest.TestCase):
    def setUp(self):
        # parse_input(), ideally for a short time window so it doesn't take too long
        f.parse_input('fluxes/sample-gamma.txt', 'e', 108, 115)
        
    def test_prepare_evt_gen(self):
        # 1) should return None
        # 2) check that it's actually doing nothing? e.g. compare dNLdE before and after?
        self.assertIsNone(f.prepare_evt_gen( list(range(108,115)) ))


class NuEmissionTest(unittest.TestCase):
    def setUp(self):
        # parse input for a short time window, to run quickly
        f.parse_input('fluxes/sample-gamma.txt', 'e', 108, 115)
        f.prepare_evt_gen(list(range(108,115)))

    def test_nu_emission(self):
        # check return values for a couple of (e_Nu, time) values
        test_values = (
                       ( 3, 109, 7.188380206236757e+51),
                       (10, 109, 1.460833123430615e+53),
                       (50, 109, 4.840627808753079e+49),
                       ( 3, 113, 6.814493577220609e+51),
                       (10, 113, 1.537299529344787e+53),
                       (50, 113, 6.431698455670024e+49),
                      )
        for (eNu, time, result) in test_values:
            # allow per mille-level changes (if internals of interpolation
            # algorithm change between different versions of scipy)
            self.assertAlmostEqual(f.nu_emission(eNu, time), result, delta=2e-3*result)


if __name__ == '__main__':
    unittest.main()
