__author__ = 'benjamin'
__copyright__ = 'copyleft'

import structures
import utils

import numpy as np

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

    # Get coordinates of barycenter and computes distances during trajectory
    t.get_barycenter()

    # Compute graphs without MSF correction
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'no_msf')
    utils.graph_bary_dist(t.bary_dist, xrange(len(t.frame)), 'no_msf')

    # Get mean position for each atom
    t.get_mean_pos()

    # Computes MSF
    t.get_msf_index()
    t.update_msf_domain()
    t.update_angles()

    # Get barycenter distance with MSF correction
    t.get_barycenter()

    # Compute graphs with MSF correction
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'msf')
    utils.graph_bary_dist(t.bary_dist, xrange(len(t.frame)), 'msf')
