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

import unittest, os
from tempfile import NamedTemporaryFile
from numpy import testing as nptest

import numpy as np

from bem import Electrodes, Mesh, Grid, Configuration, Result


class SimpletrapNewCase(unittest.TestCase):
    """simple trap triangulation and bem simulation"""
    def setUp(self):
        prefix = "examples/SimpleTrap/SimpleTrap"
        self.ele = Electrodes.from_trap(open("%s.ele" % prefix), scale=40e-6)
        self.tri = Mesh.from_electrodes(self.ele)
        self.tri.triangulate(opts="qQ")
        self.job = Configuration.select(self.tri, "RF").next()
        s, n = .1, 2*10
        self.grid = Grid(center=(0, 0, 1.5), step=(s, s, s), shape=(n, n, n))

    def assertBetween(self, a, b, c):
        self.assertGreaterEqual(a, b)
        self.assertLess(a, c)

    def test_triangulation(self):
        self.assertBetween(self.tri.points.shape[0], 60, 90)
        self.assertEqual(self.tri.points.shape[1], 3)
        self.assertBetween(self.tri.triangles.shape[0], 40, 70)
        self.assertEqual(self.tri.triangles.shape[1], 3)

    def test_select(self):
        self.assertEqual(self.job.name, "RF")
        self.assertEqual(self.tri.keys().index("RF"),
                self.job.potentials.argmax())

    def test_grid(self):
        self.assertEqual(self.grid.to_points().shape, (20**3, 3))

    def test_adapt(self):
        self.job.adapt_mesh(triangles=2e2)
        self.assertBetween(self.job.mesh.points.shape[0], 180, 220)
        self.job.adapt_mesh(triangles=1e3)
        self.assertBetween(self.job.mesh.points.shape[0], 700, 1100)

    def test_simulate(self):
        # refine twice adaptively with increasing number of triangles
        self.job.adapt_mesh(triangles=5e2)
        self.job.adapt_mesh(triangles=1e3)
        # solve for charges, get potentials and fields
        result = self.job.simulate(self.grid)

    def test_simulate_no_fmm(self):
        # refine twice adaptively with increasing number of triangles
        self.job.adapt_mesh(triangles=5e2)
        self.job.adapt_mesh(triangles=1e3)
        # solve for charges, get potentials and fields
        result = self.job.simulate(self.grid, num_lev=1)


@unittest.skip("cpy broken")
class SimpletrapCase(unittest.TestCase):
    """simple surface electrode trap with gaps"""

    def setUp(self, fil="examples/SimpleTrap/SimpleTrap.cpy"):
        mesh, mesh_data = mesh_cpy(open(fil)) 
        #os.system("gmsh %s" % (mesh.name,))

        r = 25e-6
        self.n = n = 20
        origin = [0, -r, r]
        step = [1, 2*r/n, 2*r/n]
        shape = [1, n, n]
        self.xm = grid(origin, step, shape)

        gen = simulate(mesh_data, self.xm, [("RF", True, True, True)])
        gen = list(gen)
        self.pot, self.field = gen[0][1][0][1], gen[0][1][1][1].T
        results_to_vtk("simple_trap", origin, step, shape, gen)

    def test_trap(self):
        from matplotlib import pyplot as plt
        x, n, f = self.xm, self.n, self.field
        fig, ax = plt.subplots()
        ax.contour(x[1][0], x[2][0],
                (f**2).sum(axis=0).reshape((n, n)),
                levels=np.linspace(0, 1e7, 100), cmap=plt.cm.Greys)
        #plt.show()


if __name__ == "__main__":
    unittest.main()
