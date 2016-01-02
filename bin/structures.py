__author__ = 'benjamin'

from GromacsTrajectory import XTCTrajectory
import re
import numpy as np
import mathops as mo
import copy as cp


import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class Md(object):

    def __init__(self):
        self.time = []

        self.topology = Topology()
        self.frame = []

        self.angles = []
        self.bary_dist = []

        self.mean_d1 = None
        self.mean_d2 = None

        self.msf_d1 = []
        self.msf_d2 = []


    def get_traj(self, file_input):
        traj = XTCTrajectory(file_input)
        while traj.readStep():
            f = Frame()
            f.create_frame(traj.coordinates, traj.frame_counter, traj.natoms, traj.step, traj.time)
            self.frame.append(f)

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

    def update_traj(self):
        for i in xrange(len(self.frame)):
            #print self.frame[i].xyz[-1]
            self.frame[i].update_frame(self.topology.residues, self.topology.iresidues, self.topology.atoms, self.topology.iatoms)
            #print self.frame[i].xyz[-1]

    def get_domains(self, domain_1, domain_2):
        self.topology.update_domain(domain_1)
        self.topology.update_domain(domain_2)

        self.topology.domains[0].update_eig()
        self.topology.domains[1].update_eig()

        for i in xrange(len(self.frame)):
            self.frame[i].update_domain(domain_1)
            self.frame[i].update_domain(domain_2)

            self.frame[i].domains[0].update_eig()
            self.frame[i].domains[1].update_eig()

            self.angles.append(self.frame[i].calc_angle(0,1))

    @staticmethod
    def update_domain(domain, residues, iresidues, atoms, iatoms, xyz):

        d = Domain()

        # Get indices of residues and atoms from the start to the end of the domain
        match_res = [i for i, x in enumerate(iresidues) if (x>=domain[0] and x<=domain[1])]
        match_atom = [i for i, x in enumerate(iatoms) if (x>=domain[0] and x<=domain[1])]

        # Fill residues fields : all residues between the limits of the domain
        d.residues = [residues[i] for i in match_res]
        d.iresidues = [iresidues[i] for i in match_res]
        d.nresidues = len(d.residues)

        # Fill atoms fields
        # Select all atoms between the limits of the domain
        atoms_list = [atoms[i] for i in match_atom if atoms[i]]
        # Select the indexes of the CA only
        match_ca = [i for i, x in enumerate(atoms_list) if (atoms_list[i][-2:] == 'CA')]
        d.atoms = [atoms_list[i] for i in match_ca]
        d.iatoms = [iatoms[i] for i in match_atom]
        d.natoms = len(d.atoms)


        # Select only the coordinates of atoms between the start and the end of
        # the domain
        xyz_ca = [match_atom[i] for i in match_ca]
        d.xyz = xyz[xyz_ca, :]

        # Add to topology
        return d


    def get_mean_pos(self):

        # Create a matrix that will contain the mean position for each residue for
        # both domains and initialize it at the first frame
        self.mean_d1 = self.frame[0].domains[0].xyz.copy()
        self.mean_d2 = self.frame[0].domains[1].xyz.copy()

        # Add each frame to previous matrix
        for f in xrange(1,len(self.frame)):
            self.mean_d1 = np.add(self.mean_d1, self.frame[f].domains[0].xyz.copy())
            self.mean_d2 = np.add(self.mean_d2, self.frame[f].domains[1].xyz.copy())

        # Divide by number of frames
        self.mean_d1 = self.mean_d1 / len(self.frame)
        self.mean_d2 = self.mean_d2 / len(self.frame)

    def get_msf_index(self):

        for f in xrange(len(self.frame)):

            dist_d1 = (self.mean_d1 - self.frame[f].domains[0].xyz)**2
            dist_d1 = dist_d1.sum(axis=-1)
            dist_d1 = np.sqrt(dist_d1)


            dist_d2 = (self.mean_d2 - self.frame[f].domains[1].xyz)**2
            dist_d2 = dist_d2.sum(axis=-1)
            dist_d2 = np.sqrt(dist_d2)

            if f == 0 :
                d1 = np.zeros(shape=dist_d1.shape)
                d2 = np.zeros(shape=dist_d2.shape)

            d1 = np.add(d1,dist_d1)
            d2 = np.add(d2,dist_d2)

        self.msf_d1 = cp.deepcopy(mo.box(d1))

        self.msf_d2 = cp.deepcopy(mo.box(d2))

    def update_msf_domain(self):
        for f in xrange(len(self.frame)):

            self.frame[f].domains[0].atoms = [self.frame[f].domains[0].atoms[i] for i in self.msf_d1]
            self.frame[f].domains[1].atoms = [self.frame[f].domains[1].atoms[i] for i in self.msf_d2]

            self.frame[f].domains[0].natoms = len(self.frame[f].domains[0].atoms)
            self.frame[f].domains[1].natoms = len(self.frame[f].domains[0].atoms)

            self.frame[f].domains[0].xyz = self.frame[f].domains[0].xyz[self.msf_d1,:]
            self.frame[f].domains[1].xyz = self.frame[f].domains[1].xyz[self.msf_d2,:]

    def update_angles(self):

        self.angles = []

        for i in xrange(len(self.frame)):

            self.frame[i].domains[0].eig = Eig()
            self.frame[i].domains[1].eig = Eig()

            self.frame[i].domains[0].update_eig()
            self.frame[i].domains[1].update_eig()

            self.angles.append(self.frame[i].calc_angle(0,1))

    def get_barycenter(self):

        self.bary_dist = []

        for f in xrange(len(self.frame)):
            self.frame[f].domains[0].calc_barycenter()
            self.frame[f].domains[1].calc_barycenter()

            print self.frame[f].domains[0].barycenter, self.frame[f].domains[1].barycenter

            dist = np.linalg.norm(self.frame[f].domains[0].barycenter - self.frame[f].domains[1].barycenter )
            self.bary_dist.append(dist)


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
        self.domains.append(Md.update_domain(domain, self.residues, self.iresidues, self.atoms, self.iatoms, self.xyz))


class Frame(object):

    def __init__(self):

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.domains = []

    def create_frame(self, coordinates, frame_counter, natoms, step, time):

        self.xyz = coordinates.copy()
        self.frame_counter = frame_counter
        self.step = step
        self.time = time
        self.natoms = natoms

    def update_frame(self, residues, iresidues, atoms, iatoms):

        self.residues = residues
        self.iresidues  = iresidues
        self.nresidues = len(residues)

        self.atoms = atoms
        self.iatoms = iatoms
        self.natoms = len(atoms)

    def update_domain(self, domain):
        self.domains.append(Md.update_domain(domain, self.residues, self.iresidues, self.atoms, self.iatoms, self.xyz))

    def calc_angle(self, i, j):
        #vec = mo.orientation(self.domains[i].eig.first_eig, self.domains[j].eig.first_eig)
        #return mo.angular(vec, self.domains[j].eig.first_eig)
        #return mo.angular(self.domains[i].eig.first_eig, self.domains[j].eig.first_eig)
        return mo.angular2(self.domains[i].eig.first_eig, self.domains[j].eig.first_eig)


class Domain(object):
    def __init__(self):

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0,3))

        self.eig = Eig()
        self.barycenter = []


    def update_eig(self):

        self.eig.xyz_centered = mo.center(self.xyz)
        self.eig.inertia = mo.inertia(self.eig.xyz_centered)
        self.eig.eigenvalues = mo.eigen(self.eig.inertia)[0]
        self.eig.eigenvectors = mo.eigen(self.eig.inertia)[1]

        self.eig.first_eig = self.eig.eigenvectors[:,2].tolist()

    def calc_barycenter(self):

        mw = {'ALA': 71.04, 'CYS': 103.01, 'ASP': 115.03, 'GLU': 129.04,
              'PHE': 147.07, 'GLY': 57.02, 'HIS': 137.06, 'ILE': 113.08,
              'LYS': 128.09, 'LEU': 113.08, 'MET': 131.04, 'ASN': 114.04,
              'PRO': 97.05, 'GLN': 128.06, 'ARG': 156.10, 'SER': 87.03,
              'THR': 101.05, 'VAL': 99.07, 'TRP': 186.08, 'TYR': 163.06 }

        self.barycenter = mo.barycenter(mw, self.atoms, self.xyz)

    def plane(self):
        mn = np.min(self.xyz, axis=0)
        mx = np.max(self.xyz, axis=0)

        X,Y = np.meshgrid(np.arange(mn[0], mx[0]), np.arange(mn[1], mx[1]))
        XX = X.flatten()
        YY = Y.flatten()

        # best-fit linear plane
        A = np.c_[self.xyz[:,0], self.xyz[:,1], np.ones(self.xyz.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, self.xyz[:,2]) # coefficients
        print C
        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]

        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

        # plot points and fitted surface
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(self.xyz[:,0], self.xyz[:,1], self.xyz[:,2], c='r', s=50)
        plt.xlabel('X')
        plt.ylabel('Y')
        ax.set_zlabel('Z')
        ax.axis('equal')
        ax.axis('tight')
        plt.show()


class Eig(object):
    def __init__(self):
        self.xyz_centered = np.empty(shape=(0,3))
        self.inertia = np.empty(shape=(0,3))
        self.eigenvalues = []
        self.eigenvectors = np.empty(shape=(3,3))
        self.first_eig = []


