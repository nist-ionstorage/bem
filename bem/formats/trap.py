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

from collections import OrderedDict
import logging, itertools

import numpy as np

# from enthought.tvtk.api import tvtk

class Electrodes(OrderedDict):
    """
    A set of named electrodes, each consisting of multiple faces
    (planar), each face consisting of multiple signed loops:
    {name: [[(sign, coordinates), ...more loops], ...more faces]}
    Faces are counterclockwise, negative loops are holes.
    """
    @classmethod
    def from_trap(cls, trap, scale=1.):
        """load ed electrodes in 'trap' format (inventor export)"""
        electrodes = cls()
        state = None
        for line in trap:
            line = line.rstrip()
            if not line.strip() or line.startswith("#"):
                pass
            elif line.startswith("WP"):
                pass # ignored
            elif line.startswith("S,"):
                state = "face"
                name = line.split(", ")[1]
                if name.startswith("TRAPELECTRODE_"):
                    name = name[len("TRAPELECTRODE_"):]
                if not name in electrodes:
                    electrodes[name] = []
            elif line.startswith("OUTERLOOP"):
                state = "loop"
                coords = []
                face.append((1, coords))
            elif line.startswith("INNERLOOP"):
                state = "loop"
                coords = []
                face.append((-1, coords))
            elif line.startswith("ENDOFFACE"):
                state = "face"
            else:
                #elif line.startswith(" ") or line.startswith("-") or \
                #    line.startswith("0."):
                try:
                    vector = map(float, line.split())
                    if len(vector) != 3:
                        raise ValueError("wrong length")
                except ValueError as e:
                    logging.warn("failed to parse line in trap %s: %s",
                            trap, line)
                if state == "face":
                    face = []
                    #face.normal = vector
                    electrodes[name].append(face)
                else:
                    if coords and np.allclose(coords[-1], vector):
                        logging.warn("two subsequent points are close %s, %s",
                                coords[-1], vector)
                    else:
                        coords.append(vector)

        # cleanup
        for name, faces in electrodes.iteritems():
            for face in faces[:]:
                for loop in face[:]:
                    face.remove(loop)
                    sign, coords = loop
                    coords = np.array(coords)/scale
                    if coords.shape[0] < 3:
                        logging.warn("not a loop: %s %s, removing", coords, name)
                    else:
                        face.append((sign, coords))
                if not face:
                    logging.warn("empty face: %s %s, removing", face, name)
                    faces.remove(face)

        return electrodes

    @classmethod
    def from_polygons(cls, poly):
        """load loops from 2d polygons (shapely package)"""
        obj = cls()
        for name, mp in poly:
            try:
                mp = list(mp)
            except TypeError:
                mp = [mp]
            face = []
            for pi in mp:
                ext = pi.exterior
                if not ext.is_ccw:
                    ext.coords = list(ext.coords[::-1])
                co = np.array(ext.coords[:-1])
                loop = np.zeros((co.shape[0], 3))
                loop[:, :2] = co
                face.append((1, loop))
                for ii in pi.interiors:
                    if not ii.is_ccw:
                        ii.coords = list(ii.coords[::-1])
                    co = np.array(ii.coords[:-1])
                    loop = np.zeros((co.shape[0], 3))
                    loop[:, :2] = co
                    face.append((-1, loop))
            obj[name] = [face]
        return obj

    @classmethod
    def from_system(cls, sys):
        """load loops from a planar 2d gapless electrode system
        (electrode package)"""
        # FIXME: should add ground plane where there are not electrodes
        obj = cls()
        for ele in sys:
            face = []
            for sign, coords in zip(ele.orientations(), ele.paths):
                loop = np.zeros((coords.shape[0], 3))
                loop[:, :2] = coords
                face.append((sign, loop))
            obj[ele.name] = [face]
        return obj

    def to_vtk(self, filename):
        raise NotImplemented

    def cleanup(self, tol=1e-9):
        """remove close adjacent points"""
        for name, faces in self.iteritems():
            for face in faces:
                for i, loop in enumerate(face[:]):
                    sign, coords = loop
                    coords_next = np.roll(coords, 1, axis=0)
                    l = np.square(coords-coords_next).sum(axis=1)
                    coords = coords[l > tol**2]
                    face[i] = (sign, coords)

    def extrude(self, thickness, names=None):
        """
        extrude the (xy planar) electrode surfaces by thickness into 
        into the -z direction. adds one face per pair of loop points.
        """
        if names is None:
            names = self.keys()
        dz = np.array([0., 0, thickness])
        for name in names:
            for sign, coords in itertools.chain(*self[name]):
                for ai, bi in zip(coords, np.roll(coords, 1, axis=0)):
                    coords = np.array([ai, bi, bi+dz, ai+dz])
                    if sign < 0:
                        coords = coords[::-1]
                    self[name].append([(1, coords)])

if __name__ == "__main__":
    import sys
    ele = Electrodes.from_trap(open(sys.argv[1]))
    #print ele
