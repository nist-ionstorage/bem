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

import os, re, operator
import multiprocessing
from tempfile import NamedTemporaryFile

import numpy as np

from .fastlap import TRIANGLE
from .fastlap_wrappers import sing_pot_field
from .formats.trap import Electrodes
from .formats.cpy import read_cpy
from .formats.gmsh import cpy_to_geo, read_msh, msh_to_fastlap
from .formats.vtk import results_to_vtk, read_vtk_component

default_footer = """
Field[1] = Cylinder;
Field[1].Radius = 100e-6;
Field[1].VIn = 10e-6;
Field[1].VOut = 100e-6;
Background Field = 1;
Mesh.Algorithm = 6;
"""

def mesh_cpy(fil, footer=default_footer):
    """
    takes a cpy format electrode geometry description in file object `fil`,
    converts it to gmsh geometry language, appends the gmsh footer `footer`,
    runs gmsh to mesh the geometry and returns the mesh file handle and
    the fastlap format panel data (see `simulate()` for details).
    """
    cpy = read_cpy(fil)
    geo = NamedTemporaryFile(suffix=".geo")
    cpy_to_geo(geo, cpy)
    geo.write(footer)
    geo.flush()
    mesh = NamedTemporaryFile(suffix=".msh")
    os.system("gmsh -o %s -2 %s" % (mesh.name, geo.name))
    geo.close()
    # os.system("gmsh %s" % (mesh.name,))
    points, panels, names = read_msh(mesh)
    # mesh.close()
    x, shape, groups = msh_to_fastlap(points, panels)
    return mesh, (x, shape, groups, names["name"])

def simulate(mesh, xm, selected, **opts):
    """performs the simulations for given mesh data `mesh`, observation
    points `xm` and electrode selection `selected`. other options passed
    down to `sing_pot_field()`.

    `mesh` is the tuple containing
      * the (n,4,3) array of the vertices `x`
      * the (n) array defining the shape of each panel
      * the (n) array defining which group each panel belongs to
      (indices into `names`)
      * the list of electrode names `names`
    
    `selected` is a list of tuples containing
      * the name of the electrode, and a list of booleans specifying
      what to calculate:
      * potential
      * field
      * pseudopotential (field abs squared)

    for each electrode, yields electrode name and list of tuples
    (value type name, value array at xm)
    """
    x, shape, groups, names = mesh
    xm = xm.reshape((3, -1)).T
    for name, do_pot, do_field, do_pseudo in selected:
        if not (do_pot or do_field or do_pseudo):
            continue
        potx = 1.*(groups-1==list(names).index(name))
        sing, areas, pot, field = sing_pot_field(
                x, potx, xm, shape, do_field=do_field or do_pseudo,
                **opts)
        res = []
        if do_pot:
            res.append(("potential", pot))
        if do_field:
            res.append(("field", field.T))
        if do_pseudo:
            res.append(("pseudo", (field**2).sum(axis=0)))
        yield name, res

def grid(origin, step, shape):
    """returns the rectangular uniform grid with the first point at
    `origin`, the cell diagonal `step` and a shape `shape`"""
    s = [slice(o, o+s*(h-1), 1j*h) for o, s, h in zip(origin, step, shape)]
    return np.mgrid[s]

def select(names, pots, fields, do_pseudo=True):
    """
    convert the two sequences of regexps `pots` and `fields` to flag lists
    compatible for `simulate()` given the electrode `names`.
    """
    for name in names:
        do_pot = reduce(operator.or_, (bool(re.match(p, name)) for p in pots))
        do_field = reduce(operator.or_, (bool(re.match(p, name)) for p in fields))
        yield name, do_pot, do_field, do_field and do_pseudo
