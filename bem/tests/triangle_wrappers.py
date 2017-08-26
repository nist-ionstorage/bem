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

import unittest
from numpy import testing as nptest

import numpy as np
from scipy import constants as ct

from ..triangulation import ThreeTwoTransform, Triangulation


class BasicTriangleWrapperCase(unittest.TestCase):
    def test_3to2to3(self):
        x = np.arange(9).reshape(-1, 3).astype(np.double)
        x[2, 0] = 20
        t = ThreeTwoTransform.from_points(x)
        x0, i, j = t.x0, t.i, t.j
        y = t.threed_to_twod(x)
        nptest.assert_allclose(np.linalg.norm(i), 1)
        nptest.assert_allclose(np.linalg.norm(j), 1)
        nptest.assert_allclose(y.shape, (x.shape[0], 2))
        nptest.assert_allclose(i, np.array([1, 1, 1])/3**.5)
        x1 = t.twod_to_threed(y)
        nptest.assert_allclose(x, x1)
    
    def test_3to2to3_rand(self):
        x = np.random.randn(300).reshape(-1, 3)
        x[:, 0] = 99.
        t = ThreeTwoTransform.from_points(x)
        y = t.threed_to_twod(x)
        x1 = t.twod_to_threed(y)
        nptest.assert_allclose(x, x1)

    def test_run_triangulate(self):
        surf = np.random.randn(9).reshape(-1, 3)
        tri = Triangulation.from_face([(1, surf)])
        tri.triangulate(opts="Q")
        points, triangles = tri.points, tri.triangles
        nptest.assert_allclose(surf, points)
        nptest.assert_allclose(sorted(triangles[0]), [0, 1, 2])

    @unittest.skip("can see no reason this should work")
    def test_run_triangulate_2(self):
        surf = np.random.randn(300).reshape(-1, 3)
        surf[:, 2] = 10
        tri = Triangulation.from_face([(1, surf)])
        tri.triangulate(opts="Q")
        nptest.assert_allclose(surf, tri.points)

    def test_hole(self):
        outer = np.array([(0, 0, 0), (4, 0, 0), (0, 4, 0.)])
        inner = np.array([(1, 1, 0), (2, 1, 0), (1, 2, 0.)])
        tri = Triangulation.from_face([(1, outer), (-1, inner)])
        tri.triangulate(opts="Q")
        nptest.assert_allclose(tri.points,
                np.concatenate((outer, inner)))
        self.assertNotIn((3, 4, 5), map(tuple, tri.triangles))
        nptest.assert_allclose(tri.triangles,
                [[2, 0, 5], [3, 0, 4], [5, 0, 3],
                 [5, 4, 1], [1, 4, 0], [5, 1, 2]])
