

import itertools as it
import pynomial as pn
import numpy as np


verts = []
for x,y,z in it.product([-1,1], [-1,1], [-1,1]):
    verts.append([x,y,z]);
verts = np.array(verts)/2.0;
print(verts)
np.random.seed(592765);
dt = np.random.normal(0, 0.1, verts.shape);
print(dt)
verts2 = verts+dt;

p = pn.polyhedron.Polyhedron(verts)
q = pn.polyhedron.Polyhedron(verts2)
p.merge_facets();
q.merge_facets();
# sdr = pn.polyhedron._polyhedron._SDR(p)
# sdr2 = pn.polyhedron._polyhedron._SDR(q)

sigma, rp, rq = pn.polyhedron.similarity(p, q, 0.95, 10);
print("sigma = ", sigma)
print("rp = ", rp)
print("rq = ", rq)

sigma2, rp, rq = pn.polyhedron.similarity(p, q, 0.95, 100);
print("sigma = ", sigma2)
print("rp = ", rp)
print("rq = ", rq)
