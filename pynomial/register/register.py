#
#
#
from . import _register
from . import transform
import itertools

numpy = None;
try:
    import numpy
except ImportError:
    numpy = None;

def ensure_list(data):
    if numpy and type(data) == numpy.ndarray:
        list_pts = data.tolist();
    else:
        list_pts = list(data);
    return list_pts;

#
# class icp:
#     def __init__(self, data):
#         if numpy and type(data) == numpy.ndarray:
#             list_pts = data.tolist();
#         else:
#             list_pts = list(data);
#
#         if len(list_pts) < 5:
#             msg = "Current version of icp requires at least 5 points"
#             raise(RuntimeError, msg);
#
#         assert(type(list_pts) == list and type(list_pts[0]) == list)
#         self.cpp_icp = _procrustes.icp(list_pts);
#         self.Rotation = None;
#         self.Translation = None;
#         self.rmsd = None;
#
#
#     def set_max_iterations(self, n):
#         self.cpp_icp.SetMaxIterations(n);
#
#     def set_tol(self, tol):
#         self.cpp_icp.SetMinDeltaParam(tol);
#
#     def fit(self, points, trials = 1000):
#         dim = 3;
#         # for theta in numpy.linspace(0, 2*numpy.pi, trials):
#         #     ## TODO: put this into a helper function.
#         #     q = numpy.zeros((4,));
#         #     axis = numpy.random.uniform(-1, 1, (3,));
#         #     axis = axis/numpy.sqrt(numpy.dot(axis, axis));
#         #     if dim == 2:
#         #         axis = numpy.array([0.0, 0.0, 1.0]);
#         #     q[0] = numpy.cos(theta/2.0);
#         #     q[1:] = numpy.sin(theta/2.0)*axis;
#         #     tpoints = numpy.zeros(points.shape)
#         #     for i,p in enumerate(points):
#         #         if len(p) == 2:
#         #             dim = 2;
#         #             p = numpy.array([p[0], p[1], 0])
#         #         pt = transform.quatRot(q, p);
#         #         tpoints[i, 0:dim] = pt[0:dim]
#             ##################################################
#             # do the fit.
#         icp_return = self.cpp_icp.fit(ensure_list(points));
#             # save the data if it is the smallest
#             # if not self.rmsd or icp_return.rmsd < self.rmsd:
#         self.Rotation = icp_return.R;
#         self.Translation = icp_return.t;
#         self.rmsd = icp_return.rmsd;
#             # print "theta = ", theta
#         return (self.Rotation, self.Translation, self.rmsd);

class brute_force(_register.BruteForce):
    def __init__(self, data, threshold=1e-6):
        _register.BruteForce.__init__(self, data, threshold);
        self.Rotation = None;
        self.rmsd = None;
        self.data = ensure_list(data);
        # self.cpp_brute_force = _register.BruteForce(self.data, threshold);

    def set_num_iterations(self, num):
        self.SetNumShuffles(num);

    def fit(self, points):
        self.Fit(points);
        self.Rotation = self.GetRotation();
        self.rmsd = self.GetCost();

        return (self.Rotation, self.rmsd);












#end
