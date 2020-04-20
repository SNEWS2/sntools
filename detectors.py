from math import pi
import random


# Detector Materials
water = {"molecular_weight": 18.0153, # g/mol
         "density": 1.0, # g/cm^3
         "channel_weights": {"ibd": 2, "es": 10, "o16e": 1, "o16eb": 1} # targets per molecule
}

# liquid scintillator: approximated here as CH_2
ls = {"molecular_weight": 14.0266, # g/mol
      "density": 1.0, # g/cm^3
      "channel_weights": {"ibd": 2, "es": 8, "c12e": 1, "c12eb": 1, "c12nc": 1}
}

def wbls(x):
    """Generates dictionary characterizing Water-based Liquid Scintillator.

    Input: Fraction of liquid scintillator.
    Output: Dictionary, analogous to `water` and `ls` above.
    """
    mw = x * ls["molecular_weight"] + (1-x) * water["molecular_weight"]
    d = x * ls["density"] + (1-x) * water["density"]
    cw = {}
    for channel in set(list(water["channel_weights"]) + list(ls["channel_weights"])):
        weight = x * ls["channel_weights"].get(channel, 0) \
                 + (1-x) * water["channel_weights"].get(channel, 0)
        cw[channel] = weight
    return {"molecular_weight": mw, "density": d, "channel_weights": cw}


# List of supported detector configurations
supported_detectors = ["HyperK", "SuperK", "WATCHMAN", "WATCHMAN-LS", "WATCHMAN-WbLS"]

class Detector(object):
    """A neutrino detector."""
    def __init__(self, name):
        self.name = name
        if name == "HyperK": # 2019 optimized design
            self.height = 6580.
            self.radius = 6480./2
            # 2018 Design Report: radius = 7080./2; height = 5480.
            self.material = water
        elif name == "SuperK":
            self.height = 3620.
            self.radius = 3368.15/2
            self.material = water
        elif name == "WATCHMAN": # arXiv:1502.01132
            self.height = 1280.
            self.radius = 1280./2
            self.material = water
        elif name == "WATCHMAN-LS":
            self.height = 1280.
            self.radius = 1280./2
            self.material = ls
        elif name == "WATCHMAN-WbLS":
            self.height = 1280.
            self.radius = 1280./2
            self.material = wbls(0.03) # 3% LS, 97% water
        else:
            raise ValueError("Unknown detector name: %s" % name)

        # calculate number of target molecules in detector
        volume = pi * self.radius**2 * self.height
        number_density = self.material["density"] * 6.022e23 / self.material["molecular_weight"]
        self.n_molecules = volume * number_density

    def __repr__(self):
        return "Detector('%s')" % self.name

    def __setattr__(self, attr, value):
        if attr == "n_molecules" and hasattr(self, attr):
            raise AttributeError('Number of molecules is determined by detector'
                                  ' size and material. It cannot be changed.')
        object.__setattr__(self, attr, value)

    def generate_random_vertex(self):
        while True:
            x = random.uniform(-self.radius, self.radius)
            y = random.uniform(-self.radius, self.radius)
            if x**2 + y**2 < self.radius**2:
                break
        z = random.uniform(-self.height/2, self.height/2)
        return x, y, z
