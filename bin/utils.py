__author__ = 'Benjamin Millot'
__date__ = "08/01/2016"

import argparse
import sys
import matplotlib.pyplot as plt


class Parse(object):
    """ Create a Parse class that will contain all arguments provided by user
    """
    def __init__(self):
        """Request file path to start the program

        Init:
            topology: string, mandatory, path to topology file
            trajectory: string, mandatory, path to trajectory file
            firstdomain: integer, two values, mandatory, start and end of the
                         first domain
            seconddomain: integer, two values, mandatory, start and end of the
                          second domain
            output: string, optional, path to output save file

        """
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument("-top",
                                 "--topology",
                                 help="Path to the topology",
                                 required=True,
                                 type=str)

        self.parser.add_argument("-traj",
                                 "--trajectory",
                                 help="Path to the trajectory",
                                 required=True,
                                 type=str)

        self.parser.add_argument("-fd",
                                 "--firstdomain",
                                 help="Limits of the first domain",
                                 nargs=2,
                                 required=True,
                                 type=int)

        self.parser.add_argument("-sd",
                                 "--seconddomain",
                                 help="Limits of the second domain",
                                 nargs=2,
                                 required=True,
                                 type=int)

        self.parser.add_argument("-o",
                                 "--output",
                                 help="Path to save the output",
                                 default="stdout",
                                 type=str)

        self.parser.parse_args(namespace=self)

    def commandline(self):
        """ Print comand line to either stdout or log file
        """

        out = ''

        if self.output != "stdout":
            out = " -o {0}".format(self.output)

        print "\nCommand line : module.py -top {0} -traj {1} " \
              "-fd {2} {3} -sd {4} {5}{6}\n".format(self.topology,
                                                    self.trajectory,
                                                    self.firstdomain[0],
                                                    self.firstdomain[1],
                                                    self.seconddomain[0],
                                                    self.seconddomain[1],
                                                    out)

    def redirect(self):
        """ Redirect standard output to log file
        """

        if self.output != "stdout":
            sys.stdout = open(self.output, 'w')


def graph_angles(angles, frames, name):
    """ Create the graph of the angles between the two domains during the
        trajectory

        :param angles: list of integer containing angle values
        :param frames: list of the number of frame
        :param name: string, specify if corrected (msf) or not (no_msf)
    """

    plt.plot(frames, angles, 'ro', markersize=1)
    plt.title('Angle between S and P domains during the trajectory ({0})'.
              format(name))
    plt.ylabel('Angle (in degrees)')
    plt.xlabel('Step')
    plt.savefig('../graphs/angles_{0}'.format(name))
    plt.clf()


def graph_bary_dist(bary_dist, frames, name):
    """ Create the graph of the distance between the two domains during the
        trajectory

    :param bary_dist: list of integer containing the distances
    :param frames: list of the number of frame
    :param name: string, specify if corrected (msf) or not (no_msf)
    """

    plt.plot(frames, bary_dist, 'ro', markersize=1)
    plt.title('Distance between S and P domains during the trajectory ({0})'.
              format(name))
    plt.ylabel('Distance (in nanometers)')
    plt.xlabel('Step')
    plt.savefig('../graphs/bary_dist_{0}'.format(name))
    plt.clf()
