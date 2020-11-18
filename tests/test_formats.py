import unittest

import sntools.formats as f


class TimeTest(unittest.TestCase):
    def test_get_starttime(self):
        self.assertEqual(f.get_starttime(12, 3), 12)
        self.assertEqual(f.get_starttime(0, -10), 0)
        self.assertEqual(f.get_starttime(None, -10), -10)

    def test_get_starttime_raises(self):
        with self.assertRaises(ValueError):
            f.get_starttime(10, 20)

    def test_get_endtime(self):
        self.assertEqual(f.get_endtime(-75, -10), -75)
        self.assertEqual(f.get_endtime(0, 100), 0)
        self.assertEqual(f.get_endtime(None, 100), 100)

    def test_get_endtime_raises(self):
        with self.assertRaises(ValueError):
            f.get_endtime(100, 98)


if __name__ == "__main__":
    unittest.main()
