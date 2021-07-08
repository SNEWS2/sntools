# mixing parameters from P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020, 083C01 (2020)
s12 = 0.307  # sin^2 theta_12
c12 = 1 - s12
s13 = 0.0218  # sin^2 theta_13
c13 = 1 - s13

# While exiting the supernova, neutrinos experience mass hierarchy-dependent
# flavor transitions via the MSW effect. This dictionary contains 3-tuples of
#   * original flavor at production (i.e. in input files from computer simulations),
#   * mixing probability,
#   * resulting flavor in detector.
# See the section "Treatment of Neutrino Flavour Conversion" in the documentation.
transformations = {
    "noosc": (
        ("e", 1, "e"),
        ("eb", 1, "eb"),
        ("x", 2, "x"),  # scale = 2 to include both nu_mu and nu_tau
        ("xb", 2, "xb")),
    "normal": (
        ("e", s13, "e"),  # nu_e that originated as nu_e
        ("x", c13, "e"),  # nu_e that originated as nu_x
        ("eb", c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_e
        ("xb", 1 - c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_x
        ("e", c13, "x"),  # nu_x that originated as nu_e
        ("x", 1 + s13, "x"),  # nu_x that originated as nu_x
        ("eb", 1 - c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_e
        ("xb", 1 + c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_x
    ),
    "inverted": (
        ("e", s12 * c13, "e"),  # nu_e that originated as nu_e
        ("x", 1 - s12 * c13, "e"),  # nu_e that originated as nu_x
        ("eb", s13, "eb"),  # anti-nu_e that originated as anti-nu_e
        ("xb", c13, "eb"),  # anti-nu_e that originated as anti-nu_x
        ("e", 1 - s12 * c13, "x"),  # nu_x that originated as nu_e
        ("x", 1 + s12 * c13, "x"),  # nu_x that originated as nu_x
        ("eb", c13, "xb"),  # anti-nu_x that originated as anti-nu_e
        ("xb", 1 + s13, "xb"),  # anti-nu_x that originated as anti-nu_x
    ),
}


class Transformation:
    """Class describing a transformation of neutrino fluxes"""
    def __init__(self, name):
        self.name = name
        if self.name in transformations.keys():
            self.transformation = transformations[self.name]
        else:
            raise ValueError(f"Unknown transformation name: {self.name}")

    def __repr__(self):
        return f"Transformation('{self.name}')"

    def components_producing(self, flavor):
        """Yield flux components that produce the requested flavor."""
        for (original_flv, scale, detected_flv) in self.transformation:
            if detected_flv == flavor:
                yield (original_flv, scale)
