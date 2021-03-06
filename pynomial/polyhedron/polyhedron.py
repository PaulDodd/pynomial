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

    def merge_facets(self, threshold=1e-3):
        self.MergeFacets(threshold, True);
        # try:
        #     import sklearn, sklearn.cluster
        # except ImportError as err:
        #     print("merge_facets requires sklearn package. ")
        #     print("Error occured importing sklearn package:{}".format(str(err)))
        #     raise
        # d = self.dual(origin, radius);
        # verts = d.vertices
        # clusters = sklearn.cluster.DBSCAN(eps=threshold, min_samples=1).fit(verts);
        # labels = clusters.labels_
        # assert -1 not in labels
        # # Number of clusters in labels, ignoring noise if present.
        # n_clusters_ = len(set(labels));
        # if n_clusters_ < 4:
        #     assert 0
        # _verts = np.array([[0.0,0.0,0.0] for i in range(n_clusters_)]);
        # _count = [0.0]*n_clusters_
        # for i,L in enumerate(labels):
        #     _verts[L] += verts[i];
        #     _count[L] += 1.0 ;
        # for i in range(n_clusters_):
        #     _verts[i] /= _count[i]
        # pnew = Polyhedron(_verts);
        # return pnew.dual(origin, radius)

    @property
    def vertices(self):
        return self.Vertices(False, False);

    def dihedral_angles(self, merged=False):
        angles = self.GetDihedrals()
        if merged:
            mergeMap = self.GetMergedFacesMap()
        da = [];
        for i in range(len(angles)):
            f1, f2, d = (angles[i].Face0(), angles[i].Face1(), angles[i].Angle())
            if merged and mergeMap[f1] == mergeMap[f2]: # skip faces that have been merged together
                continue;
            elif merged:
                f1 = mergeMap[f1]
                f2 = mergeMap[f2]
            da.append((f1, f2, d))
        return da

    @property
    def centroid(self):
        return self.Centroid()

    @property
    def detI(self):
        return self.GetDeterminant()

    @property
    def trI(self):
        return self.GetTrace()

    @property
    def volume(self):
        return self.Volume()

    def facet_areas(self, merged = True):
        return self.FacetAreas(False);

    @property
    def IQ(self):
        SA = sum(self.facet_areas(False))
        V = self.volume;
        return 36*np.pi*V*V/SA/SA/SA;

    @property
    def curvature(self):
        return self.Curvature();

    def rotate(self, R):
        self.Transform(R);

    def facet_normals(self, merged):
        return self.GetFacetNormals(merged)

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
    return max(A.Volume() + B.Volume() - 2.0 * AB.Volume(), 0.0);

def similarity(P, Q, tolerance=0.98, iterations=-1):
    # result = _polyhedron._SimilarityResult()
    result = _polyhedron.SimilarityMeasure(P, Q, tolerance, iterations);
    return (result.sigma, result.rotationP, result.rotationQ);
