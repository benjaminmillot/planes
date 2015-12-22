__author__ = 'benjamin'

import numpy as np

def centroid(xyz):
    return xyz - xyz.mean(axis=0)

def inertia(xyz):
    return np.transpose(xyz).dot(xyz)

def eigen(xyz):
    return np.linalg.eig(xyz)

def orientation(u,v):
    dot_product = sum(p*q for p,q in zip(u,v))
    if dot_product < 0:
        u  = [x*-1 for x in vec1]
    return u

