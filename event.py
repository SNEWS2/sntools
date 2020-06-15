class Event(object):
    """A single neutrino interaction in the detector."""

    def __init__(self, code, time=None, vertex=None, incoming=None, outgoing=None):
        self.code = code  # numeric code for interaction channel
        self.time = time  # in ms
        self.vertex = vertex or (None, None, None)  # x, y, z coordinates
        self.incoming_particles = incoming or []  # list of particles, each is a tuple containing PID, energy, direction (xyz)
        self.outgoing_particles = outgoing or []

    def __repr__(self):
        return "Event(code=%i, time=%s, vertex=%s, incoming=%s, outgoing=%s)" \
                % (self.code, self.time, self.vertex, self.incoming_particles, self.outgoing_particles)

    def __setattr__(self, name, value):
        if name in ("incoming_particles", "outgoing_particles") and hasattr(self, name):
            raise AttributeError("%s is a list. Append to it instead of overwriting it." % name)
        object.__setattr__(self, name, value)

    def nuance_string(self, i):
        """Return NUANCE-formatted representation of event for writing to output file.

        Input:
            i: number of event
        Output:
            String describing event."""

        s = "$ begin\n"
        s += "$ nuance %i\n" % self.code
        s += "$ vertex %.5f %.5f %.5f %.8f\n" % (self.vertex[0], self.vertex[1], self.vertex[2], self.time)
        for (pid, e, dirx, diry, dirz) in self.incoming_particles:
            s += "$ track %i %.5f %.5f %.5f %.5f -1\n" % (pid, e, dirx, diry, dirz)
        s += "$ info 0 0 %i\n" % i
        for (pid, e, dirx, diry, dirz) in self.outgoing_particles:
            s += "$ track %i %.5f %.5f %.5f %.5f 0\n" % (pid, e, dirx, diry, dirz)
        s += "$ end\n"
        return s
