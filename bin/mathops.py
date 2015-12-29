__author__ = 'benjamin'

import numpy as np

def norm(u):
    return np.sqrt(u[0]**2+u[1]**2+u[2]**2)


def centroid(xyz):
    return xyz - xyz.mean(axis=0)


def inertia(xyz):
    return np.transpose(xyz).dot(xyz)


def eigen(xyz):
    return np.linalg.eig(xyz)


def orientation(u,v):
    dot_product = sum(p*q for p,q in zip(u,v))
    if dot_product < 0:
        u  = [x*-1 for x in u]
    return u


def angular(u,v):
    angle = round(sum(p*q for p,q in zip(u,v)) / (norm(u) * norm(v)), 5)
    angle = np.degrees(np.arccos(angle))
    return angle


def box(dat):

    size = len(dat)
    if size % 2 == 0:
        q2_index = size/2
    else :
        q2_index = (size+1)/2

    q1_index = round(q2_index - (size - q2_index)/2)
    q3_index = round(q2_index + (size - q2_index)/2)

    if q3_index >= q1_index :
        iqr = dat[q3_index] - dat[q1_index]
    else :
        iqr = dat[q1_index] - dat[q3_index]

    threshold = dat[q3_index] + 1.5*iqr

    out = [i for i, x in enumerate(dat) if (dat[i] <= threshold)]

    return out