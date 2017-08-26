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

import itertools
import numpy as np

from ..fastlap import TRIANGLE

def flatten(fun):
    def wrapper(*a, **k):
        for i in fun(*a, **k):
            for j in i:
                yield j
    return wrapper

def to_file(fun):
    def wrapper(f, *args, **kwargs):
        for i in fun(*args, **kwargs):
            print >> f, i
    return wrapper

def msh_section(name, contents):
    yield "$%s" % name
    for line in contents:
        yield line
    yield "$End%s" % name

@to_file
@flatten
def cpy_to_msh(c):
    points, panels, names = [], [], []
    i, k, l = 0, 0, 0
    for j, (n, s) in enumerate(c.iteritems()):
        names.append(n)
        for k, (p, q) in enumerate(s):
            #print n, p, q
            points.append(p)
            panels.append(np.hstack((
                (j, k+l)*np.ones((q.shape[0], 2)),
                q+i)))
            i += p.shape[0]
        l += k
    points = np.concatenate(points)
    panels = np.concatenate(panels)
    yield msh_section("MeshFormat", ["2.2 0 8"])
    yield msh_section("Nodes", [points.shape[0]]+[
        "%i %.8f %.8f %.8f" % ((i+1,)+tuple(x)) for i, x in
        enumerate(points)])
    yield msh_section("Elements", [panels.shape[0]]+[
        "%i 2 2 %i %i %i %i %i" % ((i+1,)+tuple(t)) for i, t in
        enumerate(panels)])
    yield msh_section("PhysicalNames", [len(names)]+[
        "2 %i \"%s\"" % (i+1, n) for i, n in
        enumerate(names)])
    #yield msh_section("NodeData", [])
    #yield msh_section("ElementData", [])

@to_file
def cpy_to_geo(c):
    yield "lc1 = 1e22;" # characteristic length
    names = {}
    j, k, l = 0, 0, 0
    for n, s in c.iteritems():
        for points, panels in s:
            for pi in points:
                j += 1
                yield "Point(%i) = {%g, %g, %g, lc1};" % ((j,)+tuple(pi))
            for pi in panels:
                l += 1
                for qi in range(3):
                    k += 1
                    yield "Line(%i) = {%i, %i};" % (
                            k, pi[qi]+j-len(points),
                            pi[(qi+1)%3]+j-len(points))
                yield "Line Loop(%i) = {%i, %i, %i};" % (
                        l, k-2, k-1, k)
                yield "Plane Surface(%i) = {%i};" % (l, l)
                names.setdefault(n, []).append(l)
            #l += 1
            #yield "Compound Surface(%i) = {%s};" % (l,
            #        ", ".join([str(l-m-1) for m in
            #            range(len(panels))]))
            #names.setdefault(n, []).append(l)
    for n, pi in names.iteritems():
        yield "Physical Surface(\"%s\") = {%s};" % (n,
                ", ".join(map(str, pi)))

def read_msh(f):
    for line in f:
        assert line.startswith("$"), line
        section = line[1:].strip()
        if section == "Nodes":
            p = np.loadtxt(itertools.islice(f, int(f.next())),
                dtype=[("id", np.int), ("x", np.float),
                    ("y", np.float), ("z", np.float)])
            assert p[-1][0] == p.shape[0]
            assert f.next().startswith("$End")
        elif section == "PhysicalNames":
            n = np.loadtxt(itertools.islice(f, int(f.next())),
                dtype=[("dim", np.int), ("id", np.int),
                    ("name", np.chararray)])
            n["name"] = [ni.strip("\"") for ni in n["name"]]
            assert n[-1][1] == n.shape[0]
            assert f.next().startswith("$End")
        elif section == "Elements":
            a = np.loadtxt(itertools.islice(f, int(f.next())),
                dtype=[("id", np.int), ("typ", np.int),
                ("ntag", np.int), ("phys", np.int), ("geo", np.int),
                ("a", np.int), ("b", np.int), ("c", np.int)])
            assert a[-1][0] == a.shape[0]
            assert f.next().startswith("$End")
        elif section == "MeshFormat":
            f.next()
            assert f.next().startswith("$End")
        else:
            print section
            for l in f:
                if l.startswith("$End"): break
    return p, a, n


def msh_to_fastlap(points, panels):
    p = np.array([points[i] for i in "xyz"]).T
    g = panels["phys"]
    a = np.array([panels[i] for i in "abc"]).T
    x = np.empty((a.shape[0], 4, 3), dtype=np.double)
    x[:, :3, :] = p[a-1, :]
    shape = TRIANGLE*np.ones((x.shape[0],), dtype=np.intc)
    return x, shape, g


def grid_to_gmsh_field(origin, step, dat, fil):
    """saves the structured rectangular grid `dat` starting at `origin`
    samples spaces by `step` to the file `fil`. Format is suitable for
    a gmsh "structured field"."""
    print >> fil, "%g %g %g" % origin
    print >> fil, "%g %g %g" % step
    print >> fil, "%i %i %i" % dat.shape
    np.savetxt(fil, dat.reshape((-1, dat.shape[-1])))


def main():
    import sys
    import cpy
    for f in sys.argv[1:]:
        if True:
            c = cpy.read_cpy(open(f))
            #c = cpy.concat_cpy(c)
            cpy_to_msh(open("test.msh", "w"), c)
            cpy_to_geo(open("test.geo", "w"), c)
        else:
            c = read_msh(open(f))

if __name__ == "__main__":
    main()
