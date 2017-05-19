#
#
#

from pynomial import register
import numpy as np

dim = 2;
N = 100;
pts = np.random.uniform(-1, 1, (N,3));
center=np.mean(pts, axis=0);
pts = pts - center;
pts_transformed = np.zeros((N,3));
t = np.random.uniform(-1, 1, (3,));
t = np.array([0,0,0]);
q = np.zeros((4,));
theta = np.random.uniform(0, 2.0*np.pi);
# theta = 10.0*np.pi/180.0;
axis = np.random.uniform(-1, 1, (3,));
axis = axis/np.sqrt(np.dot(axis, axis));
if dim == 2:
    axis = np.array([0.0, 0.0, 1.0]);

q[0] = np.cos(theta/2.0);
q[1:] = np.sin(theta/2.0)*axis;
for i,p in enumerate(pts):
    pts_transformed[i] = register.transform.quatRot(q, p);
#pts_transformed = pts_transformed + t;

if dim == 2:
    pts = pts[:,0:2];
    pts_transformed = pts_transformed[:,0:2];
    print(pts.shape, pts_transformed.shape)

print("pre", pts_transformed)
# pts_transformed = pts_transformed[0:8, :]
pts_transformed = np.random.permutation(pts_transformed)[0:(N-int(0.2*N)), :]
print("post", pts_transformed)

brute = register.brute_force(pts, 1e-4);
brute.set_num_iterations(1000);
R, rmsd = brute.fit(pts_transformed);
R = np.array(R);
# tret = np.array(tret);

print( "-------- Testing Input --------")
print( "Rotation: \n", q)
print( "Translation: \n", t)

print( "-------- Brute Returned --------")
print( "Rotation: \n", R)
# print( "Translation: \n", tret)
print( "RMSD: ", rmsd)

fitted_brute = np.zeros(pts_transformed.shape);
for i,p in enumerate(pts_transformed):
    fitted_brute[i] = R.dot(p);
# fitted_brute = fitted_brute;

#
# matcher = register.icp(pts);
# matcher.set_max_iterations(10000);
# matcher.set_tol(1e-6);
# R, tret, rmsd = matcher.fit(pts_transformed);
# R = np.array(R);
# tret = np.array(tret);

# print( "-------- Testing Input --------")
# print( "Rotation: \n", q)
# print( "Translation: \n", t)
#
# print( "-------- ICP Returned --------")
# print( "Rotation: \n", R)
# print( "Translation: \n", tret)
# print( "RMSD: ", rmsd)
#
# fitted = np.zeros(pts_transformed.shape);
# for i,p in enumerate(pts_transformed):
#     fitted[i] = R.dot(p);
# fitted = fitted + tret;

import matplotlib
from matplotlib import cm
from matplotlib import rc
import pylab
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
p = ax.scatter(pts[:,0], pts[:,1], c='b', alpha = 0.5, s=100);
p = ax.scatter(pts_transformed[:,0], pts_transformed[:,1], c = 'r', alpha = 0.5, s=40);
# p = ax.scatter(fitted[:,0], fitted[:,1], c = 'k', alpha = 0.9, s=20);
p = ax.scatter(fitted_brute[:,0], fitted_brute[:,1], c = 'k', alpha = 1.0, s=20);

plt.show()
