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

from ..fastlap_wrapper import sing_pot_field
from ..formats.qui import read_qui_slow


class PlatesCase(unittest.TestCase):
    """parallel plates capacitor"""

    def setUp(self, a=1e3, n=20, d=1., u=1.):
        # edge length, points per edge, separation, potential
        self.n, self.d, self.a, self.u = n, d, a, u

        xyz = np.mgrid[0:a:1j*n, 0:a:1j*n].reshape((2, -1)).T
        self.m = m = 2*(n-1)**2*2
        self.x = np.empty((m, 4, 3), np.double)
        for i in range(n-1):
            for j in range(n-1):
                self.x[2*(i*(n-1)+j), :3, :2] = (
                        xyz[n*i+j], xyz[n*i+j+1], xyz[n*(i+1)+j])
                self.x[2*(i*(n-1)+j)+1, :3, :2] = (
                        xyz[n*i+j+1], xyz[n*(i+1)+j+1], xyz[n*(i+1)+j])
        self.x[m/2:, :3] = self.x[:m/2, 2::-1]
        self.x[:m/2, :, 2] = -d/2
        self.x[m/2:, :, 2] = d/2

        self.potx = np.r_[-(u/2)*np.ones((m/2,)), (u/2)*np.ones((m/2,))]

        self.xm = np.array([[a/2., a/2., 0]])+1e-3

        self.sing, self.areas, self.centroid, self.pot, self.field = \
                sing_pot_field(self.x, self.potx, self.xm)

    def test_values(self):
        nptest.assert_allclose(self.pot, [0], atol=5e-3)
        nptest.assert_allclose(self.field.T, [[0, 0, self.u/self.d]],
                atol=5e-3, rtol=5e-3)
        c = self.a**2/self.d/(4*np.pi) # CGS units
        q = self.sing*self.areas
        q1 = q[:self.m/2].sum()
        q2 = q[self.m/2:].sum()
        nptest.assert_allclose(c*self.u, q2, rtol=5e-3)
        nptest.assert_allclose(c*self.u, -q1, rtol=5e-3)

    def test_vary(self):
        for u in .7,:
            for n in 20, 25:
                for a in 1.2e3, 1.7e3:
                    for d in 1.5, 1.8:
                        self.setUp(a, n, d, u)
                        self.test_values()


class SphereCase(unittest.TestCase):
    """unit sphere"""

    def setUp(self, fil="examples/sphere.qui"):
        x, shape, potx = read_qui_slow(open(fil))
        self.r = r = np.linspace(1.1, 20, 200)
        self.xm = xm = r[:, None]*[[1, 1, 1]]/3**.5
        self.sing, self.areas, self.centroid, self.pot, self.field = \
                sing_pot_field(x, potx, xm, shape)

    def test_values(self):
        a = self.areas.sum()
        nptest.assert_allclose(a, 4*np.pi, rtol=1e-2)
        q = (self.sing*self.areas).sum()
        nptest.assert_allclose(q, 1, rtol=1e-2)
        nptest.assert_allclose(self.pot, 1/self.r, rtol=1e-2)
        nptest.assert_allclose(self.field.T, [[1, 1,
            1]]/self.r[:, None]**2/-(3**.5), rtol=1e-2, atol=1e-2)
        nptest.assert_allclose(self.field[0], self.field[1],
                rtol=1e-2, atol=1e-2)
        nptest.assert_allclose(self.field[1], self.field[2],
                rtol=1e-2, atol=1e-2)


if __name__ == "__main__":
    unittest.main()
