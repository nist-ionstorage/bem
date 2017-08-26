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

from ..fastlap import POINT, TRIANGLE, QUADRILAT


def read_qui(f):
    """reads a qui mesh file
    only triangles (T lines) are supported
    """
    dat = np.loadtxt(f, skiprows=1,
            usecols=(2, 5, 6, 7, 8, 9, 10, 11, 12, 13))
    m = dat.shape[0]
    pot = dat[:, 0]
    x = np.empty((m, 4, 3), dtype=np.double)
    x[:, :3, :] = dat[:, 1:].reshape((-1, 3, 3))
    shape = TRIANGLE*np.ones((m,), dtype=np.intc)
    return x, shape, pot

def read_qui_slow(f):
    """reads a qui mesh file
    all line types suported
    """
    x, shape, pot = [], [], [] 
    for line in f:
        line = line.strip().split()
        typ = line.pop(0)
        if typ == "0":
            pass # comment/title
        elif typ == "P":
            x.append(map(float, line[-3:])+[0, 0, 0]*3)
            pot.append(float(line[0]))
            shape.append(POINT)
        elif typ == "T":
            x.append(map(float, line[-9:])+[0, 0, 0])
            pot.append(float(line[0]))
            shape.append(TRIANGLE)
        elif typ == "Q":
            x.append(map(float, line[-12:]))
            pot.append(float(line[0]))
            shape.append(QUADRILAT)
    x = np.array(x, np.double).reshape((-1, 4, 3))
    shape = np.array(shape, np.intc)
    pot = np.array(pot, np.double)
    return x, shape, pot
