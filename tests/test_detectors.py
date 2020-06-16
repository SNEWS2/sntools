from __future__ import division
import unittest

import detectors as d


class DetectorTest(unittest.TestCase):
    def test_init(self):
        for name in d.supported_detectors:
            detector = d.Detector(name)
            self.assertEqual(str(detector), "Detector('%s')" % name)

    def test_invalid_detector(self):
        with self.assertRaises(ValueError):
            d.Detector("INVALID_DETECTOR_NAME")

    def test_setattr(self):
        hk = d.Detector("HyperK")
        for attr in hk.__dict__:
            with self.assertRaises(AttributeError):
                setattr(hk, attr, "")

    def test_vertex_generation(self):
        hk = d.Detector("HyperK")
        x, y, z = hk.generate_random_vertex()
        self.assertLessEqual(abs(z), hk.height / 2)
        self.assertLessEqual(x ** 2 + y ** 2, hk.radius ** 2)


if __name__ == "__main__":
    unittest.main()
