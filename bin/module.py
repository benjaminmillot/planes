__author__ = 'Benjamin Millot'
__date__ = "08/01/2016"
__version__ = "1.0"
__dependencies__ = "GromacsTrajectory, re, numpy, argparse, matplotlib"

import structures
import utils

if __name__ == '__main__':

    # Parse the arguments
    arg = utils.Parse()
    arg.redirect()
    arg.commandline()

    # Create an empty Md object
    t = structures.Md()

    # Read and store topology
    print "Currently reading the topology..."
    t.get_top(arg.topology)
    print "Reading complete !\n"

    # Read and store trajectory
    print "Currently reading the trajectory..."
    t.get_traj(arg.trajectory)
    t.update_traj()
    print "Reading complete !\n"

    print "Only carbon-alpha will be taken into account for further analysis\n"

    # Get informations from the two domains
    print "Currently locating the domain you specified..."
    print "MSF refining WILL NOT be executed"
    t.get_domains(arg.firstdomain, arg.seconddomain)
    print "Domains location complete !\n"

    print "Residues that have been conserved during the trajectory are :"
    print "First Domain"
    t.frame[-1].domains[0].print_res_domain()
    print "Second Domain"
    t.frame[-1].domains[1].print_res_domain()

    print "Angle between domains has been correctly calculated !"
    t.print_angles()
    print ""

    # Get coordinates of barycenter and computes distances during trajectory
    print "Currently computing the distances between the two domains..."
    print "MSF refining WILL NOT be executed"
    t.get_barycenter()
    print t.print_dist()
    print "Done !\n"

    # Compute graphs without MSF correction
    print "Currently computing the graph of the angle between the two domains" \
          "over the trajectory"
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'no_msf')
    print "See output file in ../graphs/angles_no_msf.png " \
          "(previous files has been rewritten)\n"

    print "Currently computing the graph of the distance between the two " \
          "domains over during the trajectory..."
    utils.graph_bary_dist(t.bary_dist, xrange(len(t.frame)), 'no_msf')
    print "Done !\n"
    print "See output file in ../graphs/bary_dist_no_msf.png " \
          "(previous files has been rewritten)\n"

    # Get mean position for each atom
    t.get_mean_pos()

    # Computes MSF
    print "Mean Square Fluctuation WILL BE ENABLED for further analysis \n"
    t.get_msf_index()
    t.update_msf_domain()

    print "Residues that have been conserved during the trajectory are :"
    print "First Domain"
    t.frame[-1].domains[0].print_res_domain()
    print "Second Domain"
    t.frame[-1].domains[1].print_res_domain()

    print "Currently calculating the angles between the two domains..."
    print "MSF is enabled"
    t.update_angles()
    t.print_angles()
    print "Done !\n"

    # Get barycenter distance with MSF correction
    print "Currently calculating the distance between the two domains..."
    print "MSF is enabled"
    t.get_barycenter()
    t.print_dist()
    print "Done !\n"

    # Compute graphs with MSF correction
    print "Currently computing the graph of the angle between the two " \
          "domains over during the trajectory..."
    utils.graph_angles(t.angles, xrange(len(t.frame)), 'msf')
    print "Done !"
    print "See output file in ../graphs/angles_msf.png " \
          "(warning: previous files, if exists, have been rewritten)\n"

    print "Currently computing the graph of the distance between the two " \
          "domains over during the trajectory..."
    utils.graph_bary_dist(t.bary_dist, xrange(len(t.frame)), 'msf')
    print "Done !"
    print "See output file in ../graphs/bary_dist_msf.png " \
          "(warning : previous files, if exists, have been rewritten)\n"