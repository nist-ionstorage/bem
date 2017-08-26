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

from ..fastlap import centroid, fastlap, TRIANGLE


class BasicFastlapCase(unittest.TestCase):
    def test_centroid(self):
        x = np.array([[[0, 0, 0.], [0, 1, 2], [1, 0, 2], [0, 0, 0]]])
        c = centroid(x, TRIANGLE*np.ones((1,), np.intc))
        nptest.assert_allclose(c, np.array([[1, 1, 4]])/3.)


if __name__ == "__main__":
    unittest.main()
