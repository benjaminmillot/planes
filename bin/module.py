__author__ = 'benjamin'
__copyright__ = 'copyleft'

from GromacsTrajectory import XTCTrajectory
import structures
import utils

import numpy as np

if __name__ == '__main__':

    # Parse the arguments
    arg = utils.Parse()

    # Read the trajectory file
    traj = XTCTrajectory(arg.trajectory)

    # Creates a Md object, with informations on topology and trajectory
    t = structures.Md()
    t.get_top(arg.topology)
    t.get_traj(traj)
    t.add_info_traj()

    # Get informations from the two domains
    t.topology.update_domain(arg.firstdomain)
    t.topology.update_domain(arg.seconddomain)

