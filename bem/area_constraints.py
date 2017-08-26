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


class Constraint(object):
    def __init__(self, inside, outside=1e30):
        self.inside, self.outside = inside, outside
    def lookup(self, x):
        # return min angle, max area
        raise NotImplemented

class Cylinder(Constraint):
    def __init__(self, start, length, radius, *a, **k):
        self.start, self.length, self.radius = (start, length, radius)
        super(Cylinder, self).__init__(*a, **k)

    def lookup(self, x):
        # http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        a = self.start-self.x
        d = np.norm(np.cross(self.length, a))/np.norm(self.length+a)
        return np.where(d < self.radius, self.inside, self.outside)
        # FIXME: cylinder has infinite length

class Box(Constraint):
    def __init__(self, start=None, center=None, diagonal=None, end=None,
            *a, **k):
        if not None in (start, end):
            center = (end + start)/2.
            diagonal = np.absolute(end - start)
        else:
            assert not None in (center, diagonal)
        self.center, self.diagonal = (center, diagonal)
        super(Box, self).__init__(*a, **k)

    def lookup(self, x):
        r = np.absolute(x - self.center) < self.diagonal
        return np.where(np.all(r, axis=-1), self.inside, self.outside)

class Sphere(Constraint):
    def __init__(self, center, radius, *a, **k):
        # ellipsoidal radii
        self.center, self.radius = (center, radius)
        super(Sphere, self).__init__(*a, **k)

    def lookup(self, x):
        r = np.absolute(x - self.center)/self.radius
        np.square(r, out=r)
        return np.where(r.sum(axis=-1) < 1, self.inside, self.outside)

