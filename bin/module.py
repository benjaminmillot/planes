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

    print t.topology.xyz
    print t.frame[1].xyz

    # np.take(matrix, list_of_index_to_keep, axis=0)   # pour ne traiter qu'une partie des matrices
    # string.index('-') # pour recuperer l'index du separateur