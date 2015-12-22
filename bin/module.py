__author__ = 'benjamin'
__copyright__ = 'copyleft'

from GromacsTrajectory import XTCTrajectory
import structures
import utils

if __name__ == '__main__':

    # Parse the arguments
    arg = utils.Parse()

    # Read the trajectory file
    traj = XTCTrajectory(arg.trajectory)

    # Creates a Md object, with informations on topology and trajectory
    t = structures.Md()
    t.get_top(arg.topology)
    t.get_traj(traj)
    t.update_traj()

    # Get informations from the two domains
    t.get_domains(arg.firstdomain, arg.seconddomain)


    print t.topology.domains[1].eig.xyz_centered
    print t.topology.domains[1].eig.inertia
    print t.topology.domains[1].eig.eigenvalues
    print t.topology.domains[1].eig.eigenvectors

    print t.frame[156].domains[0].eig.eigenvectors
    print t.frame[156].domains[0].eig.eigenvectors[:,0]


