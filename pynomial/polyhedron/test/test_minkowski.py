

import itertools as it
import pynomial as pn

import numpy as np
verts = []
for x,y,z in it.product([-1,1], [-1,1], [-1,1]):
    verts.append([x,y,z]);
verts = np.array(verts)/2.0;
print(verts)

p = pn.polyhedron.Polyhedron(verts)
p.merge_facets();
sdr = pn.polyhedron._polyhedron._SDR(p)
