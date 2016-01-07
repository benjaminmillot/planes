__author__ = 'Benjamin Millot'
__date__ = "08/01/2016"

import numpy as np


def norm(u):
    """ Compute the norm of a vector

    :param u: vector
    :return: norm of the vector
    """

    return np.sqrt(u[0]**2+u[1]**2+u[2]**2)


def center(xyz):
    """ Center a xyz array

    :param xyz: numpy array Nx3 coordinates
    :return: centered array
    """

    return xyz - xyz.mean(axis=0)


def inertia(xyz):
    """ Compute the inertia matrix from xyz coordinates

    :param xyz: centered numpy array Nx3
    :return: inertia matrix
    """

    return np.transpose(xyz).dot(xyz)


def eigen(xyz):
    """ Compute eigenvectors and eigenvalues from a set of coordinates

    :param xyz: inertia matrix - numpy array Nx3
    :return: numpy eigenvalues and eigenvectors
    """

    return np.linalg.eig(xyz)


def unit(v):
    """Get unit vector related to a vector

    :param v: entry vector
    :return: unit vector
    """

    return v / np.linalg.norm(v)


def angular(v1, v2):
    """ Compute the angle in degrees between two vector in regards to their
        orientation (same direction, opposite directions, even if one of them
        equals zero)

    :param v1: first vector
    :param v2: second vector
    :return: angle between the two vectors
    """

    u = unit(v1)
    v = unit(v2)
    angle = np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1)))

    return angle


def box(dat):
    """ Compute the Box-Whisker approach from a set of data to give the
        quartiles. Values higher than (quartile 3)+1.5*IQR are considered
        outliers.
        Used with MSF to deleted residues that move a lot during the trajectory

    :param dat: list of MSF values of each residue during trajectory
    :return: index of residues to keep in dat (under the threshold)
    """

    size = len(dat)
    mtx = np.zeros(shape=(size, 2))
    mtx[:, 0] = np.sort(dat)
    mtx[:, 1] = np.argsort(dat)

    if size % 2 == 0:
        q2_index = size/2
    else:
        q2_index = (size+1)/2

    q1_index = round(q2_index - (size - q2_index)/2)
    q3_index = round(q2_index + (size - q2_index)/2)

    iqr = mtx[q3_index, 0] - mtx[q1_index, 0]

    threshold = mtx[q3_index, 0] + 1.5*iqr

    out = [int(mtx[i, 1])
           for i, x in enumerate(mtx)
           if (mtx[i, 0] <= threshold)]
    out = sorted(out)

    return out


def barycenter(mass, residues, coord):
    """ Given coordinates, mass and residues names of a set of residues,
        compute the xyz coordinates of the barycenter

    :param mass: dictionary linking residue name with molecular weight
    :param residues: list of residues names
    :param coord: numpy array Nx3 storing residues coordinates
    :return: coordinates of barycenter in numpy list [x,y,z]
    """

    bary_coord = np.empty(shape=(0, 3))
    total_mass = 0.0
    count = 0

    for res in residues:
        if res[:3] in mass.keys():
            bary_coord = np.vstack([bary_coord, mass[res[:3]]*coord[count, :]])
            total_mass += mass[res[:3]]
        count += 1

    return np.sum(bary_coord, axis=0)/total_mass
