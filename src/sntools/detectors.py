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
    "channel_weights": {"ibd": 2, "es": 8, "ps": 2, "c12e": 1, "c12eb": 1, "c12nc": 1},
}

# liquid scintillator: LAB, average structure -> C_16.65H_27.3 (C6H5CnH2n+1 where n is 95% 9-12, 5% 13-14)
lab_snoplus = {
    "molecular_weight": 227.50,  # g/mol
    "density": 0.856,  # g/cm^3
    "channel_weights": {"ibd": 27.3, "es": 127.2, "ps": 27.3, "c12e": 16.65, "c12eb": 16.65, "c12nc": 16.65},
}
lab_juno = {
    "molecular_weight": 227.50,  # g/mol                                                                                                                                             
    "density": 0.861,  # g/cm^3                                                                                                                                                       
    "channel_weights": {"ibd": 27.3, "es": 127.2, "ps": 27.3, "c12e": 16.65, "c12eb": 16.65, "c12nc": 16.65},
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
                       "WATCHMAN", "WATCHMAN-LS", "WATCHMAN-WbLS", "JUNO",
                       "THEIA25", "THEIA100", "SNOplusAV", "SNOplusEW"]


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
        elif name == "SNOplusAV": # SNO+ inner AV only
            self.shape = "sphere"
            self.radius = 600
            self.material = lab_snoplus
        elif name == "SNOplusEW": # SNO+ external water
            self.shape = "hollowSphere"
            self.innerRadius = 605
            self.outerRadius = 900
            self.material = water
        elif name == "JUNO": # JUNO central detector
            self.shape = "sphere"
            self.radius = 1770.0
            self.material = lab_juno

            
        else:
            raise ValueError(f"Unknown detector name: {name}")

        # calculate number of target molecules in detector
        if self.shape == "box":
            volume = self.x * self.y * self.z
        elif self.shape == "cylinder":
            volume = pi * self.radius ** 2 * self.height  # assumes cylindrical detector
        elif self.shape == "sphere":
            volume = (4 * pi * self.radius ** 3) / 3
        elif self.shape == "hollowSphere":
            volume = ((4 * pi * self.outerRadius ** 3) / 3) - ((4 * pi * self.innerRadius ** 3) / 3)
        else:
            raise ValueError(f"Unknown detector shape: {self.shape}")
        number_density = self.material["density"] * 6.022e23 / self.material["molecular_weight"]
        self.n_molecules = volume * number_density

    def __repr__(self):
        return f"Detector('{self.name}')"

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
        elif self.shape == "sphere":
            while True:
                x = random.uniform(-self.radius, self.radius) 
                y = random.uniform(-self.radius, self.radius)
                z = random.uniform(-self.radius, self.radius)
                if x ** 2 + y ** 2 + z ** 2 < self.radius ** 2:
                    break
        elif self.shape == "hollowSphere":
            while True:
                x = random.uniform(-self.outerRadius, self.outerRadius)
                y = random.uniform(-self.outerRadius, self.outerRadius)
                z = random.uniform(-self.outerRadius, self.outerRadius)
                if (x**2 + y**2 + z**2 < self.outerRadius**2) and (x**2 + y**2 + z**2 > self.innerRadius**2):
                    break
        else:
            raise ValueError(f"Unknown detector shape: {self.shape}")
        return x, y, z
