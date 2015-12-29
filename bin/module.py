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

    t.frame[55].domains[0].plane()
