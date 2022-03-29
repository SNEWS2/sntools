from importlib import import_module


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
# Naming of transformations follows SNEWPy conventions.
transformations = {
    "NoTransformation": (
        ("e", 1, "e"),
        ("eb", 1, "eb"),
        ("x", 2, "x"),  # scale = 2 to include both nu_mu and nu_tau
        ("xb", 2, "xb")),
    "AdiabaticMSW_NMO": (
        ("e", s13, "e"),  # nu_e that originated as nu_e
        ("x", c13, "e"),  # nu_e that originated as nu_x
        ("eb", c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_e
        ("xb", 1 - c12 * c13, "eb"),  # anti-nu_e that originated as anti-nu_x
        ("e", c13, "x"),  # nu_x that originated as nu_e
        ("x", 1 + s13, "x"),  # nu_x that originated as nu_x
        ("eb", 1 - c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_e
        ("xb", 1 + c12 * c13, "xb"),  # anti-nu_x that originated as anti-nu_x
    ),
    "AdiabaticMSW_IMO": (
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


class SNEWPYTransformation(Transformation):
    """Adapter class to turn a SNEWPY.flavor_transformation.FlavorTransformation into an sntools.transformation.Transformation

    Some limitations currently apply:
        * the FlavorTransformation must not depend on time or energy
        * it’s not possible to set additional init parameters (e.g. non-default mixing angles)
    """

    def __init__(self, name):
        from snewpy.neutrino import MassHierarchy

        self.name = name
        kwargs = {}

        if name[-4:] == "_NMO":
            self.name = name[:-4]
            kwargs['mh'] = MassHierarchy.NORMAL
        elif name[-4:] == "_IMO":
            self.name = name[:-4]
            kwargs['mh'] = MassHierarchy.INVERTED

        xf = getattr(import_module('snewpy.flavor_transformation'), self.name)(**kwargs)

        if xf.prob_ee(0, 10) != xf.prob_ee(0, 20) or xf.prob_ee(0, 10) != xf.prob_ee(0.2, 10):
            # Probability appears to depend on time and/or energy
            raise ValueError(f"The transformation 'SNEWPY-{name}' is time- or energy-dependent. This is not yet supported.")

        transformation = (
            ("e", float(xf.prob_ee(0, 0)), "e"),  # nu_e that originated as nu_e
            ("x", float(xf.prob_ex(0, 0)), "e"),  # nu_e that originated as nu_x
            ("eb", float(xf.prob_eebar(0, 0)), "eb"),  # anti-nu_e that originated as anti-nu_e
            ("xb", float(xf.prob_exbar(0, 0)), "eb"),  # anti-nu_e that originated as anti-nu_x
            ("e", 2 * float(xf.prob_xe(0, 0)), "x"),  # nu_x that originated as nu_e
            ("x", 2 * float(xf.prob_xx(0, 0)), "x"),  # nu_x that originated as nu_x
            ("eb", 2 * float(xf.prob_xebar(0, 0)), "xb"),  # anti-nu_x that originated as anti-nu_e
            ("xb", 2 * float(xf.prob_xxbar(0, 0)), "xb"),  # anti-nu_x that originated as anti-nu_x
        )

        self.transformation = tuple(c for c in transformation if c[1] > 0)

    def __repr__(self):
        return f"<Transformation object adapted from snewpy.flavor_transformation.{self.name}>"
