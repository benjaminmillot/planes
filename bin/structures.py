__author__ = 'benjamin'

import re
import numpy as np

class Md(object):
    def __init__(self):
        self.time = []           # list of number of frames
        self.frame = []
        self.topology = Topology()

    def get_traj(self, traj):
        while traj.readStep():
            self.frame.append(Frame(traj.coordinates, traj.frame_counter, traj.natoms, traj.step, traj.time))

    def get_top(self, top):
        with open(top) as gro:
            for line in gro:
                m = re.match('^(\s){1,4}([0-9]){1,4}([A-Z]){3}', line)
                if m != None:
                    field = m.string.split()
                    atom = field[0][-3:]+field[0][:-3]
                    self.topology.update_top(atom,atom+'-'+field[1])
                    self.topology.update_coordinates(float(field[3]),float(field[4]),float(field[5]))
                self.topology.update_counter()

class Topology(object):
    def __init__(self):

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0,3))

        self.domains = []

    def update_top(self, residue, atom):

        if not self.residues or self.residues[-1]!=residue:
            self.residues.append(residue)
            self.iresidues.append(int(residue[3:]))

        if not self.atoms or self.atoms[-1]!=atom:
            self.atoms.append(atom)
            self.iatoms.append(int(atom[3:atom.index('-')]))

    def update_coordinates(self,x,y,z):

        self.xyz = np.vstack([self.xyz,[x,y,z]])

    def update_counter(self):

        self.nresidues = len(self.residues)
        self.natoms = len(self.atoms)

    def update_domain(self, domain):

        d = Domain()

        # Get index of residues and atoms from the start to the end of the domain
        match_res = [i for i, x in enumerate(self.iresidues) if (x>=domain[0] and x<=domain[1])]
        match_atom = [i for i, x in enumerate(self.iatoms) if (x>=domain[0] and x<=domain[1])]

        # Fill residues fields
        d.residues = [self.residues[i] for i in match_res]
        d.iresidues = [self.iresidues[i] for i in match_res]
        d.nresidues = len(d.residues)

        # Fill atoms fields
        d.atoms = [self.atoms[i] for i in match_atom]
        d.iatoms = [self.iatoms[i] for i in match_atom]
        d.natoms = len(d.atoms)

        # Select only the coordinates of atoms between the start and the end of
        # the domain
        d.xyz = self.xyz[match_atom, :]

        # Add to topology
        self.domains.append(d)

        # test
        print super(Topology, self)

class Frame(object):

    def __init__(self, coordinates, frame_counter, natoms, step, time):

        self.xyz = coordinates
        self.frame_counter = frame_counter
        self.natoms = natoms
        self.step = step
        self.time = time
        self.domains = []

    #def update_domain(self):


class Domain(object):
    def __init__(self):

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0,3))

