#!/usr/bin/python
# -*- coding: utf-8 -*-
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


import logging, numpy as np
from multiprocessing import Pool
from scipy.constants import epsilon_0

from electrode.gds import from_gds
from electrode.polygons import (system_to_polygons, add_gaps,
        remove_overlaps, check_validity)

from bem import Electrodes, Sphere, Mesh, Grid, Simulation, Result


def main():
    prefix = "examples/system_to_3d/rfjunction"
    scale = 40e-6
    fname = "/home/rjordens/work/nist/qc-tools/trunk/bin/rfjunction_gaps.gds"

    logging.basicConfig(level=logging.DEBUG,
        format="%(asctime)s %(name)s %(processName)s[%(process)s] "
               "%(filename)s:%(funcName)s() %(levelname)s: %(message)s",
        filename="%s.log" % prefix, filemode="w")
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter("%(processName)s[%(process)s] "
        "%(filename)s:%(funcName)s(): %(message)s"))
    logging.getLogger().addHandler(sh)

    fil = open(fname)
    s = from_gds(fil, scale)
    if False:
        p = system_to_polygons(s)
        #check_validity(p)
        #p = remove_overlaps(p)
        ## p = add_gaps(p, .1)
        ele = Electrodes.from_polygons(p)
    else:
        ele = Electrodes.from_system(s)

    # ele.cleanup()
    ele.extrude(-.1)

    mesh = Mesh.from_electrodes(ele)
    mesh.triangulate(opts="Q")
    mesh.areas_from_constraints(Sphere(center=np.array([0, 0, 1.]),
               radius=3, inside=.2, outside=20.))
    mesh.triangulate(opts="qQ", new=False)
    mesh.to_vtk(prefix)

    sim = Simulation(mesh)

    #sel = list(sim.select(potentials=[], fields=["r"]))
    sel = list(sim.select(potentials=[".*"], fields=["r.*"]))

    n, s = 2*40, .05
    grid = Grid(center=(0, 0, 1.5), step=(s, s, s), shape=(n, n, n))

    pool = Pool(4)
    c = pool.map(run_job, ((sim, grid, prefix, job) for job in sel))

    print [i[0] for i in sel]
    print np.array(c)*4*np.pi*epsilon_0*scale

    Result.view(prefix, "r")


def run_job(args):
    sim, grid, prefix, job = args
    logging.info("starting job %s", job)
    mesh = sim.adapt_mesh(job, triangles=5e3, opts="qQ")
    #mesh = sim.adapt_mesh(job, triangles=10e3, mesh=mesh, opts="q28Q")
    result = sim.simulate(job, grid, mesh, fmm_on_eval=False,
        segsize=16000)
    result.to_vtk(prefix)
    return result.collect_charges()


if __name__ == "__main__":
    main()
    import wx
    wx.PySimpleApp().MainLoop()
