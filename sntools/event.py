class Event(object):
    """A single neutrino interaction in the detector."""

    def __init__(self, code, time=None, vertex=None, incoming=None, outgoing=None):
        self.code = code  # numeric code for interaction channel
        self.time = time  # in ms
        self.vertex = vertex or (None, None, None)  # x, y, z coordinates
        self.incoming_particles = incoming or []  # list of tuples containing PID, energy, direction (x, y, z)
        self.outgoing_particles = outgoing or []

    def __repr__(self):
        return f"Event(code={self.code}, time={self.time}, vertex={self.vertex}, incoming={self.incoming_particles}, outgoing={self.outgoing_particles})"

    def __setattr__(self, name, value):
        if name in ("incoming_particles", "outgoing_particles") and hasattr(self, name):
            raise AttributeError(f"{name} is a list. Append to it instead of overwriting it.")
        object.__setattr__(self, name, value)

    def nuance_string(self, i):
        """Return NUANCE-formatted representation of event for writing to output file.

        Input:
            i: number of event
        Output:
            String describing event."""

        s = "$ begin\n"
        s += f"$ nuance {self.code}\n"
        s += f"$ vertex {self.vertex[0]:.5f} {self.vertex[1]:.5f} {self.vertex[2]:.5f} {self.time:.8f}\n"
        for (pid, e, dirx, diry, dirz) in self.incoming_particles:
            s += f"$ track {pid} {e:.5f} {dirx:.5f} {diry:.5f} {dirz:.5f} -1\n"
        s += f"$ info 0 0 {i}\n"
        for (pid, e, dirx, diry, dirz) in self.outgoing_particles:
            s += f"$ track {pid} {e:.5f} {dirx:.5f} {diry:.5f} {dirz:.5f} 0\n"
        s += "$ end\n"
        return s

    def ratpac_string(self, i, events):
        """Return RAT-PAC readable HEPEVT-style representation of event for writing to output file.

        Input:
            i: number of event
            events: list of all events
        Output:
            String describing event."""

        GeV = 0.001   # convert from MeV
        mm = 10       # convert from cm
        ns = 1000000  # convert from ms

        dt = self.time
        if i > 0:
            dt -= events[i - 1].time

        s = f"{len(self.outgoing_particles)}\n"
        for idx, (pid, e, dirx, diry, dirz) in enumerate(self.outgoing_particles):
            mass = 0.0
            if pid == 11 or pid == -11:
                mass = 0.5109907
            if pid == 22:
                mass = 0.0
            if pid == 2112:
                mass = 939.56563
            p2 = (e**2) - (mass**2)
            p = p2**0.5
            px = dirx * p
            py = diry * p
            pz = dirz * p
            if idx > 0:
                dt = 0.0
            s += f"1 {pid} 0 0 {px * GeV:.8e} {py * GeV:.8e} {pz * GeV:.8e} {mass * GeV:.8e} {dt * ns:.5e} {self.vertex[0] * mm:.5e} {self.vertex[1] * mm:.5e} {self.vertex[2] * mm:.5e}\n"
        return s
