from . import _polyhedron
import numpy as np

class Polyhedron(_polyhedron.polyhedron):

    def __init__(self, verts=None, faces=None):
        _polyhedron.polyhedron.__init__(self);
        if verts is not None and faces is not None:
            self.Set(verts, faces);
        elif verts is not None:
            self.Set(verts);
        # elif faces is not None:
        #     raise RuntimeError("Faces alone cannot define a polyhedron");

    def dual(self, origin=(0,0,0), radius=1.0):
        d = Polyhedron();
        self.Dual(d, origin, radius);
        return d;


def intersection(A, B, inside=(0,0,0)):
    ret = Polyhedron()
    _polyhedron.intersection(A, B, inside, ret);
    return ret;

def distance(A, B):
    """
    must both be centered on origin.
    for unit volume 2.0(1-V(AB))
    """
    AB = intersection(A, B);
    return A.Volume() + B.Volume() - 2.0 * AB.Volume();
