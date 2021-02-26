from __future__ import division
from math import pi
import random


# Detector Materials
water = {
    "molecular_weight": 18.0153,  # g/mol
    "density": 1.0,  # g/cm^3
    "channel_weights": {"ibd": 2, "es": 10, "o16e": 1, "o16eb": 1},  # targets per molecule
}

# liquid scintillator: approximated here as CH_2
ls = {
    "molecular_weight": 14.0266,  # g/mol
    "density": 0.86,  # g/cm^3
    "channel_weights": {"ibd": 2, "es": 8, "c12e": 1, "c12eb": 1, "c12nc": 1},
}


def wbls(x):
    """Generates dictionary characterizing Water-based Liquid Scintillator.

    Input: Fraction of liquid scintillator (0 <= x <= 1).
    Output: Dictionary, analogous to `water` and `ls` above.
    """
    if not 0 <= x <= 1:
        raise ValueError("Fraction of Liquid Scintillator must be between 0 and 1!")

    mw = x * ls["molecular_weight"] + (1 - x) * water["molecular_weight"]
    d = x * ls["density"] + (1 - x) * water["density"]
    cw = {}
    for channel in set(list(water["channel_weights"]) + list(ls["channel_weights"])):
        weight = x * ls["channel_weights"].get(channel, 0) + (1 - x) * water["channel_weights"].get(channel, 0)
        cw[channel] = weight
    return {"molecular_weight": mw, "density": d, "channel_weights": cw}


# List of supported detector configurations
supported_detectors = ["HyperK", "HyperKDR", "SuperK",
                       "WATCHMAN", "WATCHMAN-LS", "WATCHMAN-WbLS",
                       "THEIA25", "THEIA100"]


class Detector(object):
    """A neutrino detector."""

    def __init__(self, name):
        self.name = name
        if name == "HyperK":  # inner detector only, 2019 optimized design
            self.shape = "cylinder"
            self.height = 6580
            self.radius = 6480 / 2
            self.material = water
        elif name == "HyperKDR":  # 2018 Design Report (outdated, only for backwards compatibility)
            self.shape = "cylinder"
            self.height = 5480
            self.radius = 7080 / 2
            self.material = water
        elif name == "SuperK":  # inner detector only
            self.shape = "cylinder"
            self.height = 3620
            self.radius = 3368.15 / 2
            self.material = water
        elif name == "WATCHMAN":  # arXiv:1502.01132
            self.shape = "cylinder"
            self.height = 1280
            self.radius = 1280 / 2
            self.material = water
        elif name == "WATCHMAN-LS":
            self.shape = "cylinder"
            self.height = 1280
            self.radius = 1280 / 2
            self.material = ls
        elif name == "WATCHMAN-WbLS":
            self.shape = "cylinder"
            self.height = 1280
            self.radius = 1280 / 2
            self.material = wbls(0.03)  # 3% LS, 97% water
        elif name == "THEIA25":  # DOI:10.1140/epjc/s10052-020-7977-8
            self.shape = "box"
            # from dimensions in paper, substract 50cm detector wall on each side
            # estimate based on discussion with M. Wurm, G. Orebi Gann
            self.x = 2000 - 100
            self.y = 1800 - 100
            self.z = 7000 - 100
            self.material = wbls(0.10)  # 10% LS, 90% water
        elif name == "THEIA100":  # dummy values resulting in 98.2 kt volume
            self.shape = "cylinder"
            self.height = 5000
            self.radius = 5000 / 2
            self.material = wbls(0.10)  # 10% LS, 90% water
        else:
            raise ValueError("Unknown detector name: %s" % name)

        # calculate number of target molecules in detector
        if self.shape == "box":
            volume = self.x * self.y * self.z
        elif self.shape == "cylinder":
            volume = pi * self.radius ** 2 * self.height  # assumes cylindrical detector
        else:
            raise ValueError("Unknown detector shape: %s" % self.shape)
        number_density = self.material["density"] * 6.022e23 / self.material["molecular_weight"]
        self.n_molecules = volume * number_density

    def __repr__(self):
        return "Detector('%s')" % self.name

    def __setattr__(self, attr, value):
        if hasattr(self, attr):
            raise AttributeError("Detector properties cannot be changed.")
        object.__setattr__(self, attr, value)

    def generate_random_vertex(self):
        if self.shape == "box":
            x = random.uniform(0, self.x)
            y = random.uniform(0, self.y)
            z = random.uniform(0, self.z)
        elif self.shape == "cylinder":
            while True:
                x = random.uniform(-self.radius, self.radius)
                y = random.uniform(-self.radius, self.radius)
                if x ** 2 + y ** 2 < self.radius ** 2:
                    break
            z = random.uniform(-self.height / 2, self.height / 2)
        else:
            raise ValueError("Unknown detector shape: %s" % self.shape)
        return x, y, z
