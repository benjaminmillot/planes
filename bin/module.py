__author__ = 'benjamin'
__wifipass__='fk4jhr9jbg'

from GromacsTrajectory import XTCTrajectory
import structures

if __name__ == '__main__':

    file_traj = "../input/traj.xtc"
    file_top = "../input/VP1.gro"

    traj = XTCTrajectory(file_traj)

    # Creates a Md object, with informations on topology and trajectory
    t = structures.Md()
    t.get_top(file_top)
    t.get_traj(traj)

    print t.topology.xyz
    print t.frame[1].xyz

    # np.take(matrix, list_of_index_to_keep, axis=0)   # pour ne traiter qu'une partie des matrices
    # string.index('-') # pour recuperer l'index du separateur