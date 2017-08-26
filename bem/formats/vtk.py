# -*- coding: utf8 -*-
#
#   bem: triangulation and fmm/bem electrostatics tools 
#
#   Copyright (C) 2011-2012 Robert Jordens <jordens@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from enthought.tvtk.api import tvtk

def cpy_to_vtk(cpy):
    pd = tvtk.AppendPolyData()
    names = []
    for name, (points, panels) in cpy.iteritems():
        p = tvtk.PolyData(points=points, polys=panels)
        names.append(name)
        p.cell_data.scalars = len(names)*np.ones(panels.shape[0])
        pd.add_input(p)
    return names, pd

def filter_cpy(c):
    names, v = cpy_to_vtk(c)
    #v = filter(v)
    v.update()
    o = v.output
    show(v)
    p = o.points.to_array()
    a = o.polys.to_array().reshape((-1, 4))[:, 1:].ravel()
    xi = p[a, :].reshape((-1, 3, 3))
    x = np.empty((xi.shape[0], 4, 3), dtype=np.double)
    x[:, :3, :] = xi
    return x

def filter(pd):
    cpd = tvtk.CleanPolyData(input=pd.output, tolerance=1e-6)
    tf = tvtk.TriangleFilter(input=cpd.output)
    #tvtk.ClipPolyData
    sf = tvtk.LinearSubdivisionFilter(
            number_of_subdivisions=3, input=tf.output)
    #df = tvtk.DecimatePro(input=sf.output, target_reduction=.9,
    #        preserve_topology=True, splitting=True, pre_split_mesh=False,
    #        boundary_vertex_deletion=True, feature_angle=1, split_angle=1)
    df = tvtk.QuadricDecimation(input=sf.output, target_reduction=.9)
    #df = tvtk.QuadricClustering(input=sf.output,
    #        division_spacing=[3e-6, 3e-6, 1e-6])
    pdf = tvtk.SmoothPolyDataFilter(input=df.output, relaxation_factor=.1,
            number_of_iterations=50, feature_angle=1, edge_angle=1,
            boundary_smoothing=False, feature_edge_smoothing=True)
    pdf.source = tf.output
    #tvtk.AppendPolyData
    return pdf

def show(d):
    l = tvtk.LookupTable(table_range=(0, 1))
    m = tvtk.PolyDataMapper(input=d.output, scalar_visibility=True,
            scalar_mode="use_cell_data")
    p = tvtk.Property(representation="s")
    a = tvtk.Actor(mapper=m, property=p)

    ren = tvtk.Renderer(background=(.1, .2, .4))
    ren.add_actor(a)
    rw = tvtk.RenderWindow(size=(600,600))
    rw.add_renderer(ren)
    rwi = tvtk.RenderWindowInteractor(render_window=rw)
    rwi.initialize()
    rwi.start()

def mesh_to_vtk(prefix, mesh):
    names, pd = cpy_to_vtk(mesh)
    pdw = tvtk.PolyDataWriter(input=pd)
    pdw.file_name = "%s_geometry.vtk" % prefix
    pdw.write()


if __name__ == "__main__":
    import sys
    import cpy
    #show(tvtk.ConeSource(resolution=36).output)
    for f in sys.argv[1:]:
        c = cpy.read_cpy(open(f))
        c = cpy.concat_cpy(c)
        n, v = cpy_to_vtk(c)
        v.update()
        o = v.output
        #print o.points.to_array(), o.polys.to_array(), o.point_data.scalars.to_array()
        v = filter(v)
        show(v)
