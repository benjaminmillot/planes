__author__ = 'benjamin'
__copyright__ = 'copyleft'

import structures
import utils


if __name__ == '__main__':

    # Parse the arguments
    arg = utils.Parse()

    # Create an empty Md object
    t = structures.Md()

    # Read and store topology
    t.get_top(arg.topology)

    # Read and store trajectory
    t.get_traj(arg.trajectory)
    t.update_traj()

    # Get informations from the two domains
    t.get_domains(arg.firstdomain, arg.seconddomain)

    # Compute graphs without MSF correction
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'no_msf')

    # Get mean position for each atom
    t.get_mean_pos()

    # Computes msf
    t.get_msf_index()

    t.update_msf_domain()

    t.update_angles()

    # Compute graphs with MSF correction
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'msf')

    for i in t.angles:
        print i

    '''
    for i in xrange(len(t.frame)):
        print t.frame[i].domains[].xyz[-1]
    '''

    '''
    print t.topology.domains[1].eig.xyz_centered
    print t.topology.domains[1].eig.inertia
    print t.topology.domains[1].eig.eigenvalues
    print t.topology.domains[1].eig.eigenvectors
    '''

    '''
    print t.frame[156].domains[0].xyz
    print t.frame[156].domains[1].xyz

    print '---------------------------'

    print t.frame[440].domains[0].xyz
    print t.frame[440].domains[1].xyz
    '''



