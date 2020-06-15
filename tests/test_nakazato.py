import unittest

from formats import nakazato as f


class ParseInputTest(unittest.TestCase):
    def test_default_startendtime(self):
        # test extracting default start- and endtime from sample input file
        starttime, endtime, _ = f.parse_input("fluxes/intp2001.data", "e", None, None)
        self.assertEqual(starttime, -50.0)
        self.assertEqual(endtime, 20000.0)

    def test_valid_startendtime(self):
        # test extracting default start- and endtime from sample input file
        starttime, endtime, times = f.parse_input("fluxes/intp2001.data", "e", 100, 119)
        self.assertEqual(starttime, 100.0)
        self.assertEqual(endtime, 119.0)
        self.assertEqual(times, [99.0, 100.0, 105.0, 110.0, 115.0, 120.0])

    def test_invalid_startendtime(self):
        # should raise ValueError if specified start-/endtime is invalid
        with self.assertRaises(ValueError):
            f.parse_input("fluxes/intp2001.data", "e", -100, None)

        with self.assertRaises(ValueError):
            f.parse_input("fluxes/intp2001.data", "e", None, 20001)

    def test_invalid_filename(self):
        # should raise IOError (FileNotFoundError in Python 3) if file doesn't exist
        with self.assertRaises(IOError):
            f.parse_input("fluxes/NONEXISTENT_FILE.data", "e", None, None)


class PrepareEvtGenTest(unittest.TestCase):
    def setUp(self):
        # parse_input(), ideally for a short time window so it doesn't take too long
        f.parse_input("fluxes/intp2001.data", "e", 100, 119)

    def test_prepare_evt_gen(self):
        # 1) should return None
        # 2) check that it's actually doing nothing? e.g. compare dNLdE before and after?
        self.assertIsNone(f.prepare_evt_gen(list(range(100, 119))))


class NuEmissionTest(unittest.TestCase):
    def setUp(self):
        # parse input for a short time window, to run quickly
        f.parse_input("fluxes/intp2001.data", "e", 100, 119)
        f.prepare_evt_gen(list(range(100, 119)))

    def test_nu_emission(self):
        # check return values for a couple of (e_Nu, time) values
        test_values = (
            (3, 100, 8.822369704285614e52),
            (10, 100, 2.259719068605022e53),
            (50, 100, 3.972374875472509e50),
            (3, 110, 8.405614242135493e52),
            (10, 110, 2.235542709449266e53),
            (50, 110, 3.911081069511243e50),
        )
        for (eNu, time, result) in test_values:
            # allow per mille-level changes (if internals of interpolation
            # algorithm change between different versions of scipy)
            self.assertAlmostEqual(f.nu_emission(eNu, time), result, delta=2e-3 * result)


if __name__ == "__main__":
    unittest.main()
