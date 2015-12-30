__author__ = 'benjamin'

import numpy as np

def norm(u):
    return np.sqrt(u[0]**2+u[1]**2+u[2]**2)


def center(xyz):
    return xyz - xyz.mean(axis=0)


def inertia(xyz):
    return np.transpose(xyz).dot(xyz)


def eigen(xyz):
    return np.linalg.eig(xyz)


def unit(v):
    return v / np.linalg.norm(v)


def angular2(v1,v2):
    u = unit(v1)
    v = unit(v2)
    angle = np.degrees(np.arccos(np.clip(np.dot(u, v),-1,1)))
    return angle

def orientation(u, v):
    #dot_product = sum(p*q for p, q in zip(u, v))
    dot_product = np.dot(u, v)
    print dot_product
    if dot_product < 0:
        u = [x*-1 for x in u]
    return u


def angular(u,v):
    angle = sum(p*q for p,q in zip(u,v)) / (norm(u) * norm(v))
    angle = np.degrees(np.arccos(angle))
    return angle


def box(dat):

    size = len(dat)
    mtx = np.zeros(shape=(size,2))
    mtx[:,0] = np.sort(dat)
    mtx[:,1] = np.argsort(dat)

    if size % 2 == 0:
        q2_index = size/2
    else :
        q2_index = (size+1)/2

    q1_index = round(q2_index - (size - q2_index)/2)
    q3_index = round(q2_index + (size - q2_index)/2)

    iqr = mtx[q3_index,0] - mtx[q1_index,0]

    threshold = mtx[q3_index,0] + 1.5*iqr

    out = [int(mtx[i,1]) for i, x in enumerate(mtx) if (mtx[i,0] <= threshold)]
    out = sorted(out)

    return out