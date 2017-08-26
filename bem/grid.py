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

from .spherical_harmonics import SphericalHarmonics


class Grid(object):
    def __init__(self, step=None, shape=None, center=None):
        """
        returns the rectangular uniform (in each dimension) grid
        with the midpoint (center of mass) at
        `center`, the cell diagonal `step` and a shape `shape`
        """
        self.step, self.shape, self.center = step, shape, center

    def get_origin(self):
        return self.center-(np.array(self.shape)-1)/2.*self.step

    def to_slices(self):
        s = [slice(c-s*(h-1)/2., c+s*(h-1)/2., 1j*h)
                for c, s, h in zip(self.center, self.step, self.shape)]
        return s

    def to_mgrid(self):
        return np.mgrid[self.to_slices()]
    
    def to_points(self):
        x = self.to_mgrid()
        x = x.reshape(x.shape[0], -1).T
        return x

    def indices_to_coordinates(self, idx):
        idx = np.atleast_1d(idx)
        return self.center+self.step*(idx-(np.array(self.shape)-1)/2.)

    def coordinates_to_indices(self, x):
        x = np.atleast_1d(x)
        return (x-np.array(self.center))/self.step+(np.array(self.shape)-1)/2.


class UnstructuredGrid(Grid):
    """An unstructured grid supplied as a (3, any) xyz array"""
    def __init__(self, points):
        points = np.atleast_2d(points)
        assert len(points.shape) == 2
        assert points.shape[0] == 3
        self.points = points
        self.shape = (points.shape[1],)
        
    def to_mgrid(self):
        raise NotImplemented
        
    def to_points(self):
        return self.points.T

 
class SphericalGrid(Grid):
    """A (`n`, 2*`n`) spherical grid with radii given in `r` around
    list-of-points `origin`.
    The `spherical_harmonic` attribute has tools for transforming
    `result.potential` data to spherical harmonics.
    """
    def __init__(self, origins, r=1e-2, n=6):
        self.spherical_harmonics = SphericalHarmonics(n)
        self.origins = np.atleast_2d(origins)
        self.r = r
        self.theta, self.phi = self.spherical_harmonics.spherical_grid()
        xyz = np.array(self.spherical_harmonics.cartesian_grid())
        self.xyz = (self.origins.T[:, :, None, None] +
                self.r * xyz[:, None, ...])
        self.shape = self.xyz.shape[1:]
        
    def to_mgrid(self):
        return self.xyz
