import uproot

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


class EventWriter:
    def __init__(self, format, outfile):
        self.format = format

        # TODO: check file endings for consistency with format, e.g. .root for ROOT_JUNO

        if format in ('NUANCE', 'RATPAC'):
            self.outfile = open(outfile, "w")
        elif format == 'ROOT_JUNO':
            self.outfile = uproot.recreate(outfile)

    def write_preamble(self, contents: str):
        """Write preamble to output file. Not supported by all formats."""
        if self.format in ('NUANCE', 'RATPAC'):
            for line in contents.splitlines():
                self.outfile.write(f"# {line}\n")

    def write_events(self, events):
        """Write list of Event objects to output file."""
        match self.format:
            case 'NUANCE':
                self._write_nuance_events(events)
            case 'RATPAC':
                self._write_ratpac_events(events)
            case 'ROOT_JUNO':
                self._write_juno_events(events)
        
        self.outfile.close()

    def _write_nuance_events(self, events):
        """Return NUANCE-formatted representation of event for writing to output file.

        Input:
            evt: Event object
            i: number of event
        Output:
            String describing event."""

        for (i, evt) in enumerate(events):
            s = "$ begin\n"
            s += f"$ nuance {evt.code}\n"
            s += f"$ vertex {evt.vertex[0]:.5f} {evt.vertex[1]:.5f} {evt.vertex[2]:.5f} {evt.time:.8f}\n"
            for (pid, e, dirx, diry, dirz) in evt.incoming_particles:
                s += f"$ track {pid} {e:.5f} {dirx:.5f} {diry:.5f} {dirz:.5f} -1\n"
            s += f"$ info 0 0 {i}\n"
            for (pid, e, dirx, diry, dirz) in evt.outgoing_particles:
                s += f"$ track {pid} {e:.5f} {dirx:.5f} {diry:.5f} {dirz:.5f} 0\n"
            s += "$ end\n"

            self.outfile.write(s)
        
        self.outfile.write("$ stop\n")

    def _write_ratpac_events(self, events):
        """Return RAT-PAC readable HEPEVT-style representation of event for writing to output file.

        Input:
            i: number of event
            events: list of all events
        Output:
            String describing event."""

        GeV = 0.001   # convert from MeV
        mm = 10       # convert from cm
        ns = 1000000  # convert from ms

        for (i, evt) in enumerate(events):
            dt = evt.time
            if i > 0:
                dt -= events[i - 1].time

            s = f"{len(evt.outgoing_particles)}\n"
            for idx, (pid, e, dirx, diry, dirz) in enumerate(evt.outgoing_particles):
                if pid == 11 or pid == -11:
                    mass = 0.5109907
                elif pid == 2112:
                    mass = 939.56563
                elif pid == 2212:
                    mass = 938.27205
                else:
                    mass = 0.0
                p = (e**2 - mass**2)**0.5
                px = dirx * p
                py = diry * p
                pz = dirz * p
                if idx > 0:
                    dt = 0.0
                s += f"1 {pid} 0 0 {px * GeV:.8e} {py * GeV:.8e} {pz * GeV:.8e} {mass * GeV:.8e} {dt * ns:.5e} {evt.vertex[0] * mm:.5e} {evt.vertex[1] * mm:.5e} {evt.vertex[2] * mm:.5e}\n"
            
            self.outfile.write(s)

    def _write_juno_events(self, events):

        class EVENT():
            def __init__(self):
                self.nparticles = 0
                self.t = [0,0]
                self.px = [0,0]
                self.py = [0,0]
                self.pz = [0,0]
                self.nuE = 0
                self.m = [0,0]
                self.pdgid = [0,0]
                self.origPDGID = 0
                self.channel = 0

        self.outfile.mktree("SNEvents",{"nparticles": "uint64", "origPDGID":"int32", "nuE":"double", "pdgid": ("int32",(2,)),"t": ("float64",(2,)),
                                        "px": ("float64",(2,)),"py":("float64",(2,)),"pz":("float64",(2,)),"m":("float64",(2,)), "channel": "int64"})

        for orig_evt in events:
            evt = EVENT()

            for idx, (pid, e, dirx, diry, dirz) in enumerate(orig_evt.outgoing_particles):
                if pid == 11 or pid == -11:
                    mass = 0.5109907
                elif pid == 2112:
                    mass = 939.56563
                elif pid == 2212:
                    mass = 938.27205
                else:
                    mass = 0.0
                p = (e**2 - mass**2)**0.5
                evt.px[idx] = dirx * p
                evt.py[idx] = diry * p
                evt.pz[idx] = dirz * p
                evt.m[idx] = mass
                evt.pdgid[idx] = pid

            if len(orig_evt.outgoing_particles) < 2:
                # is elastic scattering, second particle is a neutrino, not visible 
                evt.px[1] = 0
                evt.py[1] = 0
                evt.pz[1] = 0
                evt.m[1] = 0
                evt.pdgid[1] = 0

            evt.nparticles = len(orig_evt.outgoing_particles)
            evt.nuE = orig_evt.incoming_particles[0][1]
            evt.t = [orig_evt.time*1e6, 0]
            evt.origPDGID = orig_evt.incoming_particles[0][0]
            evt.channel = orig_evt.code
            
            self.outfile["SNEvents"].extend({"pdgid": [evt.pdgid],"px":[evt.px],"py":[evt.py],"pz":[evt.pz],"t":[evt.t],"m":[evt.m],
                                "nuE":[evt.nuE], "nparticles":[evt.nparticles], "origPDGID":[evt.origPDGID], "channel":[evt.channel]})
