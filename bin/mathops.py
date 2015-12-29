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
    angle = sum(p*q for p,q in zip(u,v)) / (norm(u) * norm(v))
    angle = np.degrees(np.arccos(angle))
    return angle
