class Event:
    """A single neutrino interaction in the detector."""
    def __init__(self, channel = None):
        self.channel = channel # ibd, es, o16e, o16eb
        self.time = None # in ms
        self.vertex = (None, None, None) # x, y, z coordinates
        self.incoming_particles = [] # list of particles, each containing pid, energy, direction (xyz)
        self.outgoing_particles = []

    def set_time(self, t):
        """Set time of event (in ms)."""
        self.time = t

    def set_vertex(self, x, y, z):
        """Set x,y,z coordinates of event vertex (in cm)."""
        self.vertex = (x, y, z)

    def add_incoming_particle(self, p):
        """Add particle p to list of incoming particles."""
        self.incoming_particles.append(p)

    def add_outgoing_particle(self, p):
        """Add particle p to list of outgoing particles."""
        self.outgoing_particles.append(p)

    def nuance_string(self, i):
        """Return NUANCE-formatted representation of event for writing to output file.

        Input:
            i: number of event
        Output:
            String describing event."""

        s = "$ begin\n$ nuance 0\n"
        s += "$ vertex %.5f %.5f %.5f %.5f\n" % (self.vertex[0], self.vertex[1], self.vertex[2], self.time)
        for (pid, e, dirx, diry, dirz) in self.incoming_particles:
            s += "$ track %i %.5f %.5f %.5f %.5f -1\n" % (pid, e, dirx, diry, dirz)
        s += "$ info 0 0 %i\n" % i
        for (pid, e, dirx, diry, dirz) in self.outgoing_particles:
            s += "$ track %i %.5f %.5f %.5f %.5f 0\n" % (pid, e, dirx, diry, dirz)
        s += "$ end\n"
        return s
