import sys

from tvtk.api import tvtk

from bem.triangulation import Mesh
from bem.formats.stl import read_stl, stl_to_mesh, check_normals
from bem.formats.cpy import split_by_normal

s = tvtk.STLReader(file_name=sys.argv[1])
s.update()
m, d = Mesh.from_polydata(s.output)
#print m, d
#print m.points
#print m.triangles
#print m.points[m.triangles, :]

r = read_stl(open(sys.argv[1], "rb"))
check_normals(*r[:2])
r = stl_to_mesh(*r)
del r["stl_0"]
print r.keys()
r = split_by_normal(r)
m = Mesh.from_mesh(r)
m.triangulate("qQ")
m.to_vtk("test")
#print m
