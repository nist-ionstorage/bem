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
import numpy as np

def read_cpy(f, scale=1.):
    """reads cpy text file and returns a dict
    {name: [(points, triangles), ...]}
    
    for each name, a list of surfaces consisting of a points array and a
    (n,3) triangles array with indices into points
    
    * only triangles are supported
    * origin (TrapCenter or Center Point) is ignored
    * the TRAPELECTRODE_ in name is stripped
    """
    surfs = OrderedDict()
    for line in f:
        line = [field.strip() for field in line.split(",")]
        if not line:
            continue
        cmd = line.pop(0)
        if cmd == "#": # comment
            pass # comment
            # print "ignoring line", cmd, line
        elif cmd == "WP":
            point_name = line.pop(0)
            if point_name.startswith("TrapCenter"):
                origin = map(float, line) # unused
            elif point_name.startswith("Center Point"):
                origin = map(float, line) # unused
            else:
                print "ignoring line", cmd, point_name, line
        elif cmd == "S":
            points, panels = [], []
            name = line.pop(0)
            if name.startswith("TRAPELECTRODE_"):
                name = name[len("TRAPELECTRODE_"):]
        elif cmd == "V":
            points.append(map(float, line))
        elif cmd == "T":
            panels.append(map(int, line))
        elif cmd == "SEND":
            points = np.array(points, dtype=np.double)
            panels = np.array(panels, dtype=np.intc)
            points /= scale
            panels -= 1 # move to 0-index
            surfs.setdefault(name, []).append((points, panels))
        else:
            print "ignoring line", cmd, line
    return surfs

def concat_cpy(cpy):
    """concatenate surfaces (points and triangle lists) of the same name
    such that a {name: (points, triangles)} remains. in place"""
    for name, dat in cpy.iteritems():
        points, panels = [], []
        n = 0
        for p, q in dat:
            panels.append(q+n-1)
            n += p.shape[0]
            points.append(p)
        cpy[name] = np.concatenate(points), np.concatenate(panels)
    return cpy

def split_by_normal(cpy):
    """split curved faces into one face per triangle (aka split by
    normal, planarize). in place"""
    for name, faces in cpy.iteritems():
        new_faces = []
        for points, triangles in faces:
            x = points[triangles, :]
            normals = np.cross(x[:, 1]-x[:, 0], x[:, 2]-x[:, 0])
            normals /= np.sqrt(np.sum(np.square(normals), axis=1))[:, None]
            if np.allclose(normals, normals[0][None, :]):
                new_faces.append((points, triangles))
            else:
                for triangle in triangles:
                    new_faces.append((points[triangle, :],
                        np.arange(3, dtype=np.intc).reshape((1, 3))))
        cpy[name] = new_faces
    return cpy

if __name__ == "__main__":
    import sys
    for f in sys.argv[1:]:
        print read_cpy(open(f))
