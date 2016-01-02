__author__ = 'benjamin'

from GromacsTrajectory import XTCTrajectory
import re
import numpy as np
import mathops as mo


class Md(object):
    """ Creates a Md object that will contain all informations required for
        the analysis
    """

    def __init__(self):
        """ Initialize an empty Md object

        Init:
            topology: Topology object
            frame: list that contain a succession of Frame objects
            angles: list that contain the angles between the two domains
            bary_dist: list that contain the distance between the two domains
        """

        self.topology = Topology()
        self.frame = []

        self.angles = []
        self.bary_dist = []

        self.mean_d1 = None
        self.mean_d2 = None

        self.msf_d1 = []
        self.msf_d2 = []

    def get_traj(self, file_input):
        """ Read and store the trajectory

        :param file_input: binary file read by Konrad Hinsen's package
                            GromacsTrajectory
        """

        traj = XTCTrajectory(file_input)
        while traj.readStep():
            f = Frame()
            f.create_frame(traj.coordinates,
                           traj.frame_counter,
                           traj.natoms,
                           traj.step,
                           traj.time)
            self.frame.append(f)

    def get_top(self, top):
        """ Read, parse and store the topology file

        :param top: topology file (.gro)
        """

        with open(top) as gro:

            for line in gro:
                m = re.match('^(\s){1,4}([0-9]){1,4}([A-Z]){3}', line)

                if m is not None:
                    field = m.string.split()
                    atom = field[0][-3:]+field[0][:-3]
                    self.topology.update_top(atom, atom+'-'+field[1])
                    self.topology.update_coordinates(float(field[3]),
                                                     float(field[4]),
                                                     float(field[5]))
                self.topology.update_counter()

    def update_traj(self):
        """Update the trajectory according to the residues informations
        """

        for i in xrange(len(self.frame)):
            self.frame[i].update_frame(self.topology.residues,
                                       self.topology.iresidues,
                                       self.topology.atoms,
                                       self.topology.iatoms)

    def get_domains(self, domain_1, domain_2):
        """ For each frame, select the two domain and computes their eigenvalues
            and eigenvectors, then compute the angle between the two domains

        :param domain_1: Domain-type object
        :param domain_2: Domain-type object
        """

        self.topology.update_domain(domain_1)
        self.topology.update_domain(domain_2)

        self.topology.domains[0].update_eig()
        self.topology.domains[1].update_eig()

        for i in xrange(len(self.frame)):
            self.frame[i].update_domain(domain_1)
            self.frame[i].update_domain(domain_2)

            self.frame[i].domains[0].update_eig()
            self.frame[i].domains[1].update_eig()

            self.angles.append(self.frame[i].calc_angle(0, 1))

    @staticmethod
    def update_domain(domain, residues, iresidues, atoms, iatoms, xyz):
        """ Create a Domain object with required informations
            Only CA-atoms of the domains will be stored

        :param domain: list of two integers with the start and the end of the
                        domain ex [1058, 1230]
        :param residues: list containing the residues
                         ex ['ASP1059', 'GLU1060']
        :param iresidues: list containing the number of the residues
        :param atoms: list containing all the atoms
                      ['ASP1059-CA', 'ASP1059-H']
        :param iatoms: list containing the number of all the atoms
        :param xyz: numpy array Nx3 containing the coordinates of the atoms
        :return: Domain-type object
        """

        d = Domain()

        # Get indices of residues and atoms from the start to the end
        match_res = [i
                     for i, x in enumerate(iresidues)
                     if domain[0] <= x <= domain[1]]
        match_atom = [i
                      for i, x in enumerate(iatoms)
                      if domain[0] <= x <= domain[1]]

        # Fill residues fields : all residues between the limits of the domain
        d.residues = [residues[i] for i in match_res]
        d.iresidues = [iresidues[i] for i in match_res]
        d.nresidues = len(d.residues)

        # Fill atoms fields
        # Select all atoms between the limits of the domain
        atoms_list = [atoms[i] for i in match_atom if atoms[i]]
        # Select the indexes of the CA only
        match_ca = [i
                    for i, x in enumerate(atoms_list)
                    if (atoms_list[i][-2:] == 'CA')]
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
        """ Get the mean position of each residue of each domain
        """

        # Create a matrix that will contain the mean position for each residue
        # for both domains and initialize it at the first frame
        self.mean_d1 = self.frame[0].domains[0].xyz.copy()
        self.mean_d2 = self.frame[0].domains[1].xyz.copy()

        # Add each frame to previous matrix
        for f in xrange(1, len(self.frame)):
            self.mean_d1 = np.add(self.mean_d1,
                                  self.frame[f].domains[0].xyz.copy())
            self.mean_d2 = np.add(self.mean_d2,
                                  self.frame[f].domains[1].xyz.copy())

        # Divide by number of frames
        self.mean_d1 /= len(self.frame)
        self.mean_d2 /= len(self.frame)

    def get_msf_index(self):
        """ Get indices of residues to keep inside each domain of each frame
            Residues to keep are determined by Box-Whisker approach
        """

        for f in xrange(len(self.frame)):

            dist_d1 = (self.mean_d1 - self.frame[f].domains[0].xyz)**2
            dist_d1 = np.sum(dist_d1, axis=1)
            dist_d1 = np.sqrt(dist_d1)

            dist_d2 = (self.mean_d2 - self.frame[f].domains[1].xyz)**2
            dist_d2 = np.sum(dist_d2, axis=1)
            dist_d2 = np.sqrt(dist_d2)

            if f == 0:
                d1 = np.zeros(shape=dist_d1.shape)
                d2 = np.zeros(shape=dist_d2.shape)

            d1 = np.add(d1, dist_d1)
            d2 = np.add(d2, dist_d2)

        self.msf_d1 = mo.box(d1)
        self.msf_d2 = mo.box(d2)

    def update_msf_domain(self):
        """ Use stored indices to update each domain of each frames
            In each domain, residue with high mobility will be deleted in terms
            of coordinates and atom names
        """

        for f in xrange(len(self.frame)):

            self.frame[f].domains[0].atoms = [self.frame[f].domains[0].atoms[i]
                                              for i in self.msf_d1]
            self.frame[f].domains[1].atoms = [self.frame[f].domains[1].atoms[i]
                                              for i in self.msf_d2]

            self.frame[f].domains[0].natoms = \
                len(self.frame[f].domains[0].atoms)
            self.frame[f].domains[1].natoms = \
                len(self.frame[f].domains[0].atoms)

            self.frame[f].domains[0].xyz = \
                self.frame[f].domains[0].xyz[self.msf_d1, :]
            self.frame[f].domains[1].xyz = \
                self.frame[f].domains[1].xyz[self.msf_d2, :]

    def update_angles(self):
        """ Compute angles between the two domain in each frame of the
            trajectory
        """

        self.angles = []

        for i in xrange(len(self.frame)):

            self.frame[i].domains[0].eig = Eig()
            self.frame[i].domains[1].eig = Eig()

            self.frame[i].domains[0].update_eig()
            self.frame[i].domains[1].update_eig()

            self.angles.append(self.frame[i].calc_angle(0, 1))

    def get_barycenter(self):
        """ Compute distance between the two domain in each frame of the
            trajectory
        """

        self.bary_dist = []

        for f in xrange(len(self.frame)):
            self.frame[f].domains[0].calc_barycenter()
            self.frame[f].domains[1].calc_barycenter()

            dist = np.linalg.norm(self.frame[f].domains[0].barycenter -
                                  self.frame[f].domains[1].barycenter)

            self.bary_dist.append(dist)


class Topology(object):
    """ Create a topology object
        Contain informations from the topology file.
        Necessary to link with the trajectory.
    """

    def __init__(self):
        """ Initialize an empty Topology object

        Init:
            residue: list of all residue ex : ['ASP1059', 'GLU1060']
            nresidue: integer, number of residue in the protein
            iresidue: list of all number of residues ex : [1059, 1060]

            atoms: list of all atoms ex : ['ASP1059-CA', 'ASP1059-H']
            natoms: integer, number of atoms in the protein
            iatoms: list of all number of atoms
       """

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0, 3))

        self.domains = []

    def update_top(self, residue, atom):
        """ Update topology with incoming informations

        :param residue: string containing the residue name ex 'ASP1059'
        :param atom: string containing the atom name ex 'ASP1059-CA'
        """

        if not self.residues or self.residues[-1] != residue:
            self.residues.append(residue)
            self.iresidues.append(int(residue[3:]))

        if not self.atoms or self.atoms[-1] != atom:
            self.atoms.append(atom)
            self.iatoms.append(int(atom[3:atom.index('-')]))

    def update_coordinates(self, x, y, z):
        """ Increment the matrix of coordinates by a row

        :param x: float, x coordinate
        :param y: float, y coordinate
        :param z: float, z coordinate
        """

        self.xyz = np.vstack([self.xyz, [x, y, z]])

    def update_counter(self):
        """ Update the total number of residue and atom after reading
        """

        self.nresidues = len(self.residues)
        self.natoms = len(self.atoms)

    def update_domain(self, domain):
        """ Add a domain to the topology

        :param domain: Domain-type object
        """

        self.domains.append(Md.update_domain(domain,
                                             self.residues,
                                             self.iresidues,
                                             self.atoms,
                                             self.iatoms,
                                             self.xyz))


class Frame(object):
    """ Create a Frame object
        Contain the informations at each step of the trajectory
    """

    def __init__(self):
        """ Initialize an empty Frame object

        Init:
            residue: list of all residue ex : ['ASP1059', 'GLU1060']
            nresidue: integer, number of residue in the protein
            iresidue: list of all number of residues ex : [1059, 1060]

            atoms: list of all atoms ex : ['ASP1059-CA', 'ASP1059-H']
            natoms: integer, number of atoms in the protein
            iatoms: list of all number of atoms

            xyz: numpy array Nx3 containing the coordinates of atoms
            frame_counter: integer, number of the frame
            time: float, time equivalent to the frame

            domains: list containing two Domain-type objects
       """

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0, 3))
        self.frame_counter = 0
        self.time = 0.0
        self.step = []

        self.domains = []

    def create_frame(self, coordinates, frame_counter, natoms, step, time):
        """ Fill informations in the current frame

        :param coordinates: numpy array Nx3 containing the coordinates of atoms
        :param frame_counter: integer, number of the frame
        :param natoms: integer, number of atoms of the protein
        :param time: float, time equivalent to the frame
        """

        self.xyz = coordinates.copy()
        self.frame_counter = frame_counter
        self.step = step
        self.time = time
        self.natoms = natoms

    def update_frame(self, residues, iresidues, atoms, iatoms):
        """ Update frame with incoming informations

            :param residues: list of all residue ex : ['ASP1059', 'GLU1060']
            :param iresidues: list of all number of residues ex : [1059, 1060]

            :param atoms: list of all atoms ex : ['ASP1059-CA', 'ASP1059-H']
            :param iatoms: list of all number of atoms
        """

        self.residues = residues
        self.iresidues = iresidues
        self.nresidues = len(residues)

        self.atoms = atoms
        self.iatoms = iatoms
        self.natoms = len(atoms)

    def update_domain(self, domain):
        """ Add a domain to the current frame

        :param domain: Domain-type object
        """

        self.domains.append(Md.update_domain(domain,
                                             self.residues,
                                             self.iresidues,
                                             self.atoms,
                                             self.iatoms,
                                             self.xyz))

    def calc_angle(self, i, j):
        """ Compute the angle between the two domains of the frame

        :param i: third eigenvector of the first domain
        :param j: third eigenvector of the second domain
        """
        return mo.angular(self.domains[i].eig.third_eig,
                          self.domains[j].eig.third_eig)


class Domain(object):
    """ Create a Domain object with informations relative to a specific domain
    """

    def __init__(self):
        """ Initialize an empty Domain object

        Init:
            residue: list of all residue in the domain
                     ex : ['ASP1059', 'GLU1060']
            nresidue: integer, number of residue in the domain
            iresidue: list of  number of residues in the domain
                      ex : [1059, 1060]

            atoms: list of atoms in the domains, restricted to CA
                   ex : ['ASP1059-CA', 'ASP1059-H']
            natoms: integer, number of atoms in the domain
            iatoms: list of the number of atoms

            xyz: numpy array Nx3 containing the coordinates of atoms

            eig: Eig-type object
            barycenter: list of the x, y, z coordinates of the barycenter of
                        the domain ex [3.53, 7.34, 4.98]
       """

        self.residues = []
        self.nresidues = 0
        self.iresidues = []

        self.atoms = []
        self.natoms = 0
        self.iatoms = []

        self.xyz = np.empty(shape=(0, 3))

        self.eig = Eig()
        self.barycenter = []

    def update_eig(self):
        """ Update the Eig object with a centered matrix in order to get the
            inertia matrix. Then compute the eigenvectors and store the third.
        """

        self.eig.xyz_centered = mo.center(self.xyz)
        self.eig.inertia = mo.inertia(self.eig.xyz_centered)
        self.eig.eigenvalues = mo.eigen(self.eig.inertia)[0]
        self.eig.eigenvectors = mo.eigen(self.eig.inertia)[1]

        self.eig.third_eig = self.eig.eigenvectors[:, 2].tolist()

    def calc_barycenter(self):
        """ Compute the distance between the two barycenter of the two domains
        """

        mw = {'ALA': 71.04, 'CYS': 103.01, 'ASP': 115.03, 'GLU': 129.04,
              'PHE': 147.07, 'GLY': 57.02, 'HIS': 137.06, 'ILE': 113.08,
              'LYS': 128.09, 'LEU': 113.08, 'MET': 131.04, 'ASN': 114.04,
              'PRO': 97.05, 'GLN': 128.06, 'ARG': 156.10, 'SER': 87.03,
              'THR': 101.05, 'VAL': 99.07, 'TRP': 186.08, 'TYR': 163.06}

        self.barycenter = mo.barycenter(mw, self.atoms, self.xyz)


class Eig(object):
    """ Create an Eig object that will contain all informations relative to the
        computation of eigenvectors
    """
    def __init__(self):
        """ Create an empty Eig object

        Init:
            xyz_centered: numpy array Nx3 containing centered positions
            inertia: numpy array Nx3 containing the inertia matrix
            eigenvalues: list of the three eigenvalues
                         ex [9.746, 13.3732, 6.74832]
            eigenvectors: numpy array 3x3 containing the three eigenvectors
            third_eig: list of the coordinates of the least eigenvector
                        ex [0.5675, -0.6532, 1.3234]
        """
        self.xyz_centered = np.empty(shape=(0, 3))
        self.inertia = np.empty(shape=(0, 3))
        self.eigenvalues = []
        self.eigenvectors = np.empty(shape=(3, 3))
        self.third_eig = []