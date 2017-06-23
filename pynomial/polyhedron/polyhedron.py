from . import _polyhedron
import numpy as np
import copy

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

    def merge_facets(self, origin=(0,0,0), radius=1.0, threshold=1e-3):
        try:
            import sklearn, sklearn.cluster
        except ImportError as err:
            print("merge_facets requires sklearn package. ")
            print("Error occured importing sklearn package:{}".format(str(err)))
            raise
        d = self.dual(origin, radius);
        verts = d.vertices
        clusters = sklearn.cluster.DBSCAN(eps=threshold, min_samples=1).fit(verts);
        labels = clusters.labels_
        assert -1 not in labels
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels));
        if n_clusters_ < 4:
            assert 0
        _verts = np.array([[0.0,0.0,0.0] for i in range(n_clusters_)]);
        _count = [0.0]*n_clusters_
        for i,L in enumerate(labels):
            _verts[L] += verts[i];
            _count[L] += 1.0 ;
        for i in range(n_clusters_):
            _verts[i] /= _count[i]
        pnew = Polyhedron(_verts);
        return pnew.dual(origin, radius)

    @property
    def vertices(self):
        return self.Vertices(False, False);

    @property
    def dihedral_angles(self):
        angles = self.GetDihedrals()
        return [(angles[i].Face0(), angles[i].Face1(), angles[i].Angle()) for i in range(len(angles))]

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
