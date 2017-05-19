#
# some of the following functions are taken from hpmc.
#
from . import _register
import numpy as np

# Multiply two quaternions
# Apply quaternion multiplication per http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
# (requires numpy)
# \param q1 quaternion
# \param q2 quaternion
# \returns q1*q2
def quatMult(q1, q2):
    s = q1[0]
    v = q1[1:]
    t = q2[0]
    w = q2[1:]
    q = np.empty((4,), dtype=float)
    q[0] = s*t - np.dot(v, w)
    q[1:] = s*w + t*v + np.cross(v,w)
    return q

def quatConj(q):
    return np.array([ q[0], -q[1], -q[2], -q[3]] );

# Rotate a vector by a unit quaternion
# Quaternion rotation per http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
# (requires numpy)
# \param q rotation quaternion
# \param v 3d vector to be rotated
# \returns q*v*q^{-1}
def quatRot(q, v):
    v = np.asarray(v)
    q = np.asarray(q)
    # assume q is a unit quaternion
    w = q[0]
    r = q[1:]
    vnew = np.empty((3,), dtype=v.dtype)
    vnew = v + 2*np.cross(r, np.cross(r,v) + w*v)
    return vnew

def translate(t, v):
    t = np.asarray(t);
    v = np.asarray(v);
    return v + t;

def ensure_list(data):
    if np and type(data) == np.ndarray:
        list_pts = data.tolist();
    else:
        list_pts = list(data);
    return list_pts;

def quat_from_matrix(rotation):
    return np.array(_register.make_quaternion_from_rotation_matrix(ensure_list(rotation)));

#
# TODO: RotationMatrixToQuaternion(R)
# TODO: QuaternionToRotationMatrix(q)
# TODO: RotationQuaternion(v, theta)
#
