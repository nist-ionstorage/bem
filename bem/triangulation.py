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

import logging, itertools
from collections import OrderedDict

import numpy as np
try:
    from tvtk.api import tvtk
except ImportError:
    import warnings
    warnings.warn("not tvtk found")


from .pytriangle import triangulate


class ThreeTwoTransform(object):
    """
    transforms planar face coordinates between 3D cartesian and 2D local
    cartesian face coordinates
    """
    x0 = None
    i = None
    j = None

    @classmethod
    def from_points(cls, x):
        """
        finds the 3D origin x[0] and
        the two 3D unit vectors i, j defining the 2D coordinate system.

        the plane is defined by the first three vertices if possible, else
        the next vertices are used.
        """
        assert x.shape[1] == 3
        x0 = x[0]
        i = x[1] - x0
        i /= np.linalg.norm(i)
        for xi in x[2:]:
            j = xi - x0
            ij = np.dot(i, j)
            j -= i*ij
            jj = np.linalg.norm(j)
            if jj > 1e-3*ij:
                # j is significantly different from i
                j /= jj
                break
        obj = cls()
        obj.x0, obj.i, obj.j = x0, i, j
        # flip axes if polygon in twod has negative area (is cw)
        # otherwise triangulate() will flip it
        # this is not overly expensive
        if polygon_area(obj.threed_to_twod(x)) < 0:
            obj.i, obj.j = obj.j, obj.i
        logging.debug("local coordinate system is x0=%s, i=%s, j=%s",
                obj.x0, obj.i, obj.j)
        return obj

    def threed_to_twod(self, x):
        """
        transforms planar surface in 3D defined by (n, 3) verices `x` to 
        2D coordinates.

        returns the (n, 2) coordinates y

        whether the vertices actually lie in that plane is ignored
        (they are projected onto the plane).
        """
        assert x.shape[1] == 3
        y = np.dot(x - self.x0, np.array([self.i, self.j]).T)
        return y

    def twod_to_threed(self, y):
        """
        undo threed_to_twod()
        """
        assert y.shape[1] == 2
        x = self.x0 + self.i * y[:, 0, None] + self.j * y[:, 1, None]
        return x


def point_inside_polygon(p, poly):
    """
    determine if a point (x, y) is inside a given polygon or not
    Polygon is a list of (x, y) pairs.

    http://www.ariel.com.au/a/python-point-int-poly.html
    """
    x, y = p
    inside = False
    n = len(poly)
    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y) and y <= max(p1y, p2y) and x <= max(p1x, p2x):
            if p1y != p2y:
                xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x or x <= xinters:
                inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def find_point_inside_polygon(y):
    """
    returns a point within the given polygon
    TODO: use better algorithm,
    http://www.exaflop.org/docs/cgafaq/cga2.html
    """
    # here we take a point just left or right of the midpoint between
    # the first two vertices. since inner loops are rare and this only
    # needs to be done once, this is not overly expensive.
    h = (y[1] + y[0]) / 2
    hp = (y[1] - y[0])
    hp = np.array([hp[1], -hp[0]])
    hp = hp / (hp**2).sum()**.5 * (h**2).sum()**.5 * 1e-6
    if not point_inside_polygon(h + hp, y):
        hp *= -1
    return h + hp


def polygon_area(polygon):
    """
    signed area of a polygon in 2D (cw: negative, ccw: positive)
    """
    polygon = np.atleast_2d(polygon)
    x1, y1 = polygon.T
    x2, y2 = np.roll(polygon.T, 1, axis=-1)
    area = ((x2-x1)*(y2+y1)).sum()
    return area


class Triangulation(object):
    """
    Triangulation of a planar face.
    
    Can consist of several loops (boundaries), inner loops (holes). See
    from_face()
    Can also be specified by points and triangles only (no loops). See
    from_mesh()
    """
    coords = None
    points = None
    triangles = None

    @classmethod
    def from_face(cls, loops):
        """
        create a face triangulation from a list of signed loops: [(sign,
        loop), ...]. If sign>0, coords define an outer loop else an inner
        loop (a hole).
        """
        points3, points2, pointmarkers = [], [], []
        segments, segmentmarkers = [], []
        holes = []
        npoints = 0
        coords = None
        for i, (sig, x) in enumerate(loops):
            if coords is None:
                coords = ThreeTwoTransform.from_points(x)
            n = x.shape[0]
            si = npoints+np.arange(n, dtype=np.intc)
            npoints += n
            # markers start at 2 (0 and 1 reserved for inside and boundary)
            mi = (i+2)*np.ones(n, dtype=np.intc)
            pointmarkers.append(mi)
            y = coords.threed_to_twod(x)
            if polygon_area(y) < 0:
                x, y = x[::-1], y[::-1]
            points3.append(x)
            points2.append(y)
            segs = np.array((np.roll(si, 1), si))
            segments.append(np.ascontiguousarray(segs.T))
            segmentmarkers.append(mi)
            if sig < 0: # mark holes
                holes.append(find_point_inside_polygon(y))

        points3 = np.concatenate(points3)
        points2 = np.concatenate(points2)

        points3_ = coords.twod_to_threed(points2)
        if not np.allclose(points3, points3_):
            raise ValueError("not planar: %s, %s" % (points3, points3_))

        holes = np.array(holes, dtype=np.double).reshape(-1, 2)

        obj = cls()
        obj.coords = coords
        obj._args = dict(
                points=points2,
                pointmarkers=np.concatenate(pointmarkers),
                segments=np.concatenate(segments),
                segmentmarkers=np.concatenate(segmentmarkers),
                holes=holes)
        return obj

    @classmethod
    def from_mesh(cls, points_triangles):
        """
        load a planar face from a (points, triangles) pair
        """
        points, triangles = points_triangles
        obj = cls()
        obj.coords = ThreeTwoTransform.from_points(points)
        obj.triangles = triangles
        obj.points = points
        obj._args = dict(
                points=obj.coords.threed_to_twod(points),
                triangles=triangles)
        return obj

    def triangulate(self, opts, new=False):
        """
        (re) triangulate this face using options in `opts`. creates a
        new triangulation if `new` is true, else modifies `self`.
        Returns the new triangulation or self.
        """
        # logging.debug("calling triangulate()")
        ret = triangulate(opts=opts, **self._args)
        # logging.debug("done with triangulate()")
        if new:
            obj = Triangulation()
        else:
            obj = self
        # since triangle() recycles arguments, not necessary
        # obj._args.update(ret)
        obj._args = ret
        obj.coords = self.coords
        obj.points = self.coords.twod_to_threed(ret["points"])
        obj.triangles = ret["triangles"]
        return obj

    def areas_from_constraints(self, *constraints):
        """
        set area constraints such that the max area of the triangles
        created in future refinements conforms to the minimum of the
        constraints. constraints must not be empty.
        """
        areas = np.nanmin([
            c.lookup(self.points) for c in constraints], axis=0)
        # take minimum over triangle vertices
        areas = areas.take(self.triangles).min(axis=-1)
        self.set_max_areas(areas)

    def set_max_areas(self, max_areas):
        self._args["triangleareas"] = max_areas

    def plot(self, ax, axes=(0, 1), *args, **kwargs):
        """
        plot 2D projection of mesh (use given axes indices, project
        along remaining axis)
        """
        ax.triplot(self.points[:, axes[0]], self.points[:, axes[1]],
                self.triangles.copy(), *args, **kwargs)


class Mesh(OrderedDict):
    """
    A Mesh is a mapping of electrode names to list of faces
    """
    points = None
    triangles = None
    groups = None

    def check(self):
        """
        verifies that points, triangles are superficially consistent
        call gather() before this.
        """
        n = self.points.shape[0]
        assert np.all(np.unique(self.triangles) == np.arange(n)), \
                "unused points"
        m = self.triangles.shape[0]
        assert self.groups.shape[0] == m, (self.groups.shape, m)
        assert self.groups.min() == 0, self.groups.min()
        assert self.groups.max() == len(self.keys())-1, (self.groups.max(),
                len(self.keys()))
        logging.debug("superficially consistent")

    @classmethod
    def from_parts(cls, data, method):
        """
        loads a mesh description {name: [data, ...], ...}
        by calling method(data) for each face of each electrode
        """
        mesh = cls()
        for name, faces in data.items():
            mesh[name] = [method(face) for face in faces]
        return mesh

    @classmethod
    def from_electrodes(cls, electrodes):
        """
        creates a Triangulation() for each face and groups them by
        electrode name.
        electrodes should be like {name: [face, ...], ...}
        face should be like [(sign, loop), ...]
        """
        return cls.from_parts(electrodes, Triangulation.from_face)

    @classmethod
    def from_mesh(cls, mesh):
        """
        reads a free-format mesh {name: [(points, triangles), ...], ...}
        """
        return cls.from_parts(mesh, Triangulation.from_mesh)

    def triangulate(self, opts="qQ", new=False):
        """
        (re)-triangulate, create a new Mesh() if new==True,
        else modify `self`
        """
        if new:
            obj = Mesh()
        else:
            obj = self
        for name, faces in self.items():
            obj[name] = [face.triangulate(opts, new) for face in faces]
        obj.gather()
        return obj

    def gather(self):
        """
        concatenate points and triangles from the individual faces and
        generates group markers (indices into the electrode
        names/self.keys())
        """
        points, triangles, groups = [], [], []
        n = 0
        for i, (name, faces) in enumerate(self.items()):
            for face in faces:
                points.append(face.points)
                triangles.append(face.triangles+n)
                n += face.points.shape[0]
                groups.append(i*np.ones((face.triangles.shape[0]),
                    dtype=np.intc))
        self.points, self.triangles, self.groups = map(np.concatenate,
                (points, triangles, groups))
        logging.info("%i points, %i triangles from %i electrodes",
                self.points.shape[0], self.triangles.shape[0],
                len(self.items()))
        self.check()
        return self.points, self.triangles, self.groups

    def fastlap_points(self):
        """returns a point array suitable for fastlap, fastlap likes
        points clockwise, triangle and vtk like them counterclockwise"""
        x = np.empty((self.triangles.shape[0], 4, 3), dtype=np.double)
        x[:, :3, :] = self.points[self.triangles[:, ::-1], :]
        return x
   
    def areas_from_constraints(self, constraints):
        """set max triangle areas for subsequent triangulations"""
        for face in itertools.chain(*self.values()):
            face.areas_from_constraints(constraints)

    def set_max_areas(self, max_areas):
        n = 0
        # TODO: could use np.split(max_areas, np.cumsum()) here
        for face in itertools.chain(*self.values()):
            ni = face.triangles.shape[0]
            face.set_max_areas(max_areas[n:n+ni])
            n += ni

    def to_polydata(self, **kwargs):
        """
        return a vtk PolyData object for this mesh, add named cell_data arrays
        from **kwargs (need to have same length as self.triangles)
        """
        pd = tvtk.PolyData(points=self.points, polys=self.triangles)
        e = tvtk.IntArray(name="electrode_index")
        e.from_array(self.groups)
        pd.cell_data.add_array(e)
        for name, data in kwargs.items():
            if data is None:
                continue
            c = tvtk.DoubleArray(name=name)
            assert data.shape[0] == self.triangles.shape[0]
            c.from_array(data)
            pd.cell_data.add_array(c)
        return pd

    def to_vtk(self, prefix, **kwargs):
        """saves mesh as prefix_mesh.vtk"""
        pd = self.to_polydata(**kwargs)
        pdw = tvtk.PolyDataWriter(input=pd)
        pdw.file_name = "%s_mesh.vtk" % prefix
        pdw.write()
        logging.debug("written mesh to %s, polydata %s",
                pdw.file_name, kwargs.keys())
        return pd

    @classmethod
    def from_polydata(cls, pd):
        """loads mesh from polydata in pd"""
        obj = cls()
        obj.points = pd.points.to_array()
        # drop the repeated first point
        obj.triangles = pd.polys.to_array().reshape(-1, 4)[:, 1:]
        datasets = {}
        for i in range(pd.cell_data.number_of_arrays):
            name = pd.cell_data.get_array_name(i)
            da = pd.cell_data.get_array(i)
            data = da.to_array()
            if name == "electrode_index":
                obj.groups = data
            else:
                datasets[name] = data
        return obj, datasets

    @classmethod
    def from_vtk(cls, prefix):
        """loads mesh from vtk polydata"""
        pdr = tvtk.PolyDataReader()
        pdr.file_name = "%s_mesh.vtk" % prefix
        pdr.update()
        return cls.from_polydata(pdr.output)

    def plot(self, ax, *args, **kwargs):
        """plot mesh"""
        for name, faces in self.items():
            for face in faces:
                face.plot(ax, *args, **kwargs)
