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

from bem.formats.trap import Electrodes

class SystemTo3DGaps(unittest.TestCase):
    @unittest.skip("needs gds file")
    def test_add_3d(self):
        from electrode.gds import from_gds
        fname = "/home/rjordens/work/nist/qc-tools/trunk/bin/rfjunction_gaps.gds"
        fil = open(fname)
        s = from_gds(fil, 40e-6)
        e = Electrodes.from_system(s)
        e.extrude(2.)
