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

from ..pytriangle import triangulate

class BasicTriangleCase(unittest.TestCase):
    def test_basic(self):
        """http://www.cs.cmu.edu/~quake/triangle.delaunay.html
        spiral.node"""
        x = np.array([ 0    ,   0, -0.416,   0.909, -1.35 ,   0.436,
             -1.64 ,  -0.549, -1.31 ,  -1.51, -0.532,  -2.17,
              0.454,  -2.41, 1.45 ,  -2.21, 2.29 ,  -1.66,
              2.88 ,  -0.838, 3.16 ,   0.131, 3.12 ,   1.14,
              2.77 ,   2.08, 2.16 ,   2.89, 1.36 ,   3.49]).reshape(-1, 2)
        ret = triangulate(points=x, opts="zQ")
        nptest.assert_allclose(ret["pointmarkers"], 
        np.array([0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
        nptest.assert_allclose(ret["points"], 
        np.array([[ 0.   ,  0.   ], [-0.416,  0.909], [-1.35 ,  0.436],
               [-1.64 , -0.549], [-1.31 , -1.51 ], [-0.532, -2.17 ],
               [ 0.454, -2.41 ], [ 1.45 , -2.21 ], [ 2.29 , -1.66 ],
               [ 2.88 , -0.838], [ 3.16 ,  0.131], [ 3.12 ,  1.14 ],
               [ 2.77 ,  2.08 ], [ 2.16 ,  2.89 ], [ 1.36 ,  3.49 ]]))
        nptest.assert_allclose(ret["triangles"], 
        np.array([[ 0,  4,  5], [ 0,  5,  6], [ 2,  3,  0],
               [ 1,  0, 12], [ 0,  1,  2], [ 3,  4,  0],
               [ 2,  1, 14], [ 6,  7,  0], [ 8,  9,  0],
               [10,  0,  9], [10, 11,  0], [ 1, 13, 14],
               [13,  1, 12], [12,  0, 11], [ 7,  8,  0]]))

    def test_tricall(self):
        """triangle/tricall.c"""
        p = np.array([.0, .0, 1., .0, 1., 10., .0, 10.]).reshape(-1, 2)
        pa = np.array([.0, 1., 11., 10]).reshape(-1, 1)
        pm = np.array([0, 2, 0, 0], dtype=np.intc)
        r = np.array([.5, 5., 7., .1]).reshape(-1, 1)
        ret = triangulate(points=p, pointattributes=pa, pointmarkers=pm,
                regions=r, opts="pczAevnQ")
        nptest.assert_allclose(ret["edgemarkers"], np.array([1, 1, 0, 1, 1]))
        nptest.assert_allclose(ret["edges"],
            np.array([[3, 0], [0, 1], [1, 3], [1, 2], [2, 3]]))
        nptest.assert_allclose(ret["neighbors"],
            np.array([[-1,  1, -1], [-1,  0, -1]]))
        nptest.assert_allclose(ret["pointattributes"],
            np.array([[  0.], [  1.], [ 11.], [ 10.]]))
        nptest.assert_allclose(ret["pointmarkers"], np.array([1, 2, 1, 1]))
        nptest.assert_allclose(ret["points"],
            np.array([[  0.,   0.], [  1.,   0.], [  1.,  10.], [  0.,  10.]]))
        nptest.assert_allclose(ret["segmentmarkers"], np.array([1, 1, 1, 1]))
        nptest.assert_allclose(ret["segments"],
            np.array([[1, 0], [2, 1], [3, 2], [0, 3]]))
        nptest.assert_allclose(ret["triangleattributes"],
            np.array([[ 7.], [ 7.]]))
        nptest.assert_allclose(ret["triangles"],
            np.array([[3, 0, 1], [1, 2, 3]]))
        nptest.assert_allclose(ret["voroniedges"],
            np.array([[ 0, -1], [ 0, -1], [ 0,  1], [ 1, -1], [ 1, -1]]))
        nptest.assert_allclose(ret["voroninorms"],
            np.array([[-10.,   0.], [  0.,  -1.], [  0.,   0.],
                [ 10.,   0.], [  0.,   1.]]))
        nptest.assert_allclose(ret["voronipointattributes"],
            np.array([[ 5.5], [ 5.5]]))
        nptest.assert_allclose(ret["voronipoints"],
            np.array([[ 0.5,  5. ], [ 0.5,  5. ]]))
        for k in "edges edgemarkers voroniedges neighbors voroninorms "\
                 "voronipoints voronipointattributes".split():
            del ret[k]
        ta = np.array([3., 1.])
        ret = triangulate(triangleareas=ta, opts="QprazBP", **ret)
        nptest.assert_allclose(ret["pointattributes"], 
        np.array([[  0.  ], [  1.  ], [ 11.  ], [ 10.  ], [  5.  ],
               [  9.25], [  6.  ], [  8.5 ], [  3.5 ], [  4.25],
               [  7.5 ], [  6.75], [  2.5 ], [  1.75]]))
        nptest.assert_allclose(ret["points"], 
        np.array([[  0.  ,   0.  ], [  1.  ,   0.  ], [  1.  ,  10.  ],
               [  0.  ,  10.  ], [  0.  ,   5.  ], [  0.5 ,   8.75],
               [  1.  ,   5.  ], [  1.  ,   7.5 ], [  1.  ,   2.5 ],
               [  0.5 ,   3.75], [  0.  ,   7.5 ], [  0.5 ,   6.25],
               [  0.  ,   2.5 ], [  0.5 ,   1.25]]))
        nptest.assert_allclose(ret["triangleattributes"], 
        np.array([[ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.],
               [ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.], [ 7.]]))
        nptest.assert_allclose(ret["triangles"], 
        np.array([[ 1, 13,  0], [ 4, 11, 10], [ 2,  5,  7], [ 6,  4,  9],
               [ 8, 12, 13], [ 6, 11,  4], [12,  0, 13], [ 2,  3,  5],
               [ 7, 10, 11], [ 8,  6,  9], [ 4, 12,  9], [ 3, 10,  5],
               [ 7,  5, 10], [ 6,  7, 11], [ 8,  9, 12], [ 1,  8, 13]]))


if __name__ == "__main__":
    unittest.main()
