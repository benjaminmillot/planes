import argparse
import matplotlib.pyplot as plt


class Parse(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument("-top", "--topology", help="Path to the topology", required=True, type=str)
        self.parser.add_argument("-traj", "--trajectory", help="Path to the trajectory", required=True, type=str)

        self.parser.add_argument("-fd", "--firstdomain", help="Limits of the first domain", nargs=2, required=True, type=int)
        self.parser.add_argument("-sd", "--seconddomain", help="Limits of the second domain", nargs=2, required=True, type=int)

        self.parser.add_argument("-o", "--output", help="Path to save the output", default="stdout", type=str)

        self.parser.parse_args(namespace=self)


def graph_angles(angles, frames, name):

    plt.plot(frames, angles, 'ro', markersize=1)
    plt.title('Angle between S and P domains during the trajectory ({0})'.format(name))
    plt.ylabel('Angle (in degrees)')
    plt.xlabel('Step')
    plt.savefig('../graphs/angles_{0}'.format(name))
    plt.clf()

def graph_bary_dist(bary_dist, frames, name):

    plt.plot(frames, bary_dist, 'ro', markersize=1)
    plt.title('Distance between S and P domains during the trajectory ({0})'.format(name))
    plt.ylabel('Distance (in Angstrom)')
    plt.xlabel('Step')
    plt.savefig('../graphs/bary_dist_{0}'.format(name))
    plt.clf()
