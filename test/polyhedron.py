import pynomial as pn
import shapes
import numpy as np


shape = [ s  for s in shapes.SHAPES if s.ShortName == "P04"]
poly = pn.polyhedron.Polyhedron(verts=shape[0].vertices)
# poly2 = pn.polyhedron.Polyhedron(verts=shape[1].vertices)
poly.writePOS("icosahedron.pos")
dpoly = poly.dual(origin=(0.1, 0.1, 0.1))
# poly2.writePOS("octahedron.pos")
# inter = pn.polyhedron.intersection(poly, poly2);
dpoly.writePOS("dod_skew.pos")

# print(dir(poly))
# del poly

print("end test!")
