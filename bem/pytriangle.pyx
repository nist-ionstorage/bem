# -*- coding: utf8 -*-
#
#   bem: triangulation and fmm/bem electrostatics tools 
#   pytriangle - python bindings for triangle
#   http://www.cs.cmu.edu/~quake/triangle.html
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

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

import sys
import cython

import numpy as np
cimport numpy as np

np.import_array()

# NOTE: this code assumes REAL=double

# np.int_t is a long (vis python int), but ints (C) are 32 bits
ctypedef int intc_t

cdef extern from "triangle.h":
    cdef struct triangulateio:
        double *pointlist
        double *pointattributelist
        int *pointmarkerlist
        int numberofpoints
        int numberofpointattributes

        int *trianglelist
        double *triangleattributelist
        double *trianglearealist
        int *neighborlist
        int numberoftriangles
        int numberofcorners
        int numberoftriangleattributes

        int *segmentlist
        int *segmentmarkerlist
        int numberofsegments

        double *holelist
        int numberofholes

        double *regionlist
        int numberofregions

        int *edgelist
        int *edgemarkerlist
        double *normlist
        int numberofedges

    cdef void c_triangulate "triangulate"(
            char *, triangulateio *, triangulateio *,
            triangulateio *) nogil

    cdef void trifree(void *memptr) nogil


# numpy array wrapper for raw arrays that frees its underlying memory
# on numpy array deallocation
# modified from https://gist.github.com/1249305

cdef class ArrayWrapper:
    cdef void* data
    cdef set_data(self, void* data):
        self.data = data
    def __dealloc__(self):
        free(<void*>self.data)

cdef np.ndarray wrap(void* data, int dt, int dim, np.npy_intp* shape):
    cdef np.ndarray a = np.PyArray_SimpleNewFromData(dim, shape, dt, data)
    aw = ArrayWrapper()
    aw.set_data(data)
    a.base = <PyObject*> aw
    Py_INCREF(aw)
    return a

cdef PyObject *py_triunsuitable = NULL


@cython.boundscheck(False)
@cython.wraparound(False)
def triangulate(
        np.ndarray[np.double_t, ndim=2, mode="c"] points not None,
        np.ndarray[np.double_t, ndim=2, mode="c"] pointattributes=None,
        np.ndarray[intc_t, ndim=1, mode="c"] pointmarkers=None,
        np.ndarray[intc_t, ndim=2, mode="c"] triangles=None,
        np.ndarray[np.double_t, ndim=2, mode="c"] triangleattributes=None,
        np.ndarray[np.double_t, ndim=1, mode="c"] triangleareas=None,
        np.ndarray[intc_t, ndim=2, mode="c"] segments=None,
        np.ndarray[intc_t, ndim=1, mode="c"] segmentmarkers=None,
        np.ndarray[np.double_t, ndim=2, mode="c"] holes=None,
        np.ndarray[np.double_t, ndim=2, mode="c"] regions=None,
        object triunsuitable=None,
        opts="",
        ):
    """
    Parameters
    -----------

    Returns
    --------

    Other Parameters
    ----------------

    Notes
    -----

    References
    ----------

    Examples
    --------

    """
    cdef triangulateio tin, tout, tvorout
    cdef np.ndarray ndarray
    cdef char* copts
    global py_triunsuitable

    # grab the parameters and fill the three triangulateio structures
    # accordingly. also enforce necessary flags in opts
    
    opts += "z"
    tin.pointlist = &points[0, 0]
    tin.numberofpoints = points.shape[0]
    if pointattributes is not None:
        assert pointattributes.shape[0] == points.shape[0]
        tin.pointattributelist = &pointattributes[0, 0]
        tin.numberofpointattributes = pointattributes.shape[1]
    else:
        tin.numberofpointattributes = 0
    if pointmarkers is not None:
        assert pointmarkers.shape[0] == points.shape[0]
        tin.pointmarkerlist = &pointmarkers[0]
    else:
        tin.pointmarkerlist = NULL
    if triangles is not None:
        opts += "r"
        tin.trianglelist = &triangles[0, 0]
        tin.numberoftriangles = triangles.shape[0]
        tin.numberofcorners = triangles.shape[1]
        if triangleattributes is not None:
            assert triangleattributes.shape[0] == triangles.shape[0]
            tin.triangleattributelist = &triangleattributes[0, 0]
            tin.numberoftriangleattributes = triangleattributes.shape[1]
        else:
            tin.numberoftriangleattributes = 0
    else:
        tin.numberoftriangles = 0
        assert "r" not in opts
    if triangleareas is not None:
        opts += "a"
        assert triangleareas.shape[0] == triangles.shape[0]
        tin.trianglearealist = &triangleareas[0]
    else:
        #assert "a" not in opts # FIXME a0.01 is fine
        pass
    if segments is not None:
        opts += "p"
        tin.segmentlist = &segments[0, 0]
        tin.numberofsegments = segments.shape[0]
        if segmentmarkers is not None:
            assert segmentmarkers.shape[0] == segments.shape[0]
            tin.segmentmarkerlist = &segmentmarkers[0]
    else:
        tin.numberofsegments = 0
        # assert "p" not in opts
    if holes is not None:
        tin.holelist = &holes[0, 0]
        tin.numberofholes = holes.shape[0]
    else:
        tin.numberofholes = 0
    if regions is not None:
        assert regions.shape[0] == points.shape[0]
        tin.regionlist = &regions[0, 0]
        tin.numberofregions = regions.shape[1]
    else:
        tin.numberofregions = 0
 
    # let triangulate() allocate memory. we will have to free() them later
    tout.pointlist = NULL
    tout.pointattributelist = NULL
    tout.pointmarkerlist = NULL
    tout.trianglelist = NULL
    tout.triangleattributelist = NULL
    tout.neighborlist = NULL
    tout.segmentlist = NULL
    tout.segmentmarkerlist = NULL
    tout.edgelist = NULL
    tout.edgemarkerlist = NULL
    tvorout.pointlist = NULL
    tvorout.pointattributelist = NULL
    tvorout.edgelist = NULL
    tvorout.normlist = NULL

    if not triunsuitable is None:
        py_triunsuitable = <PyObject*> triunsuitable
        opts += "u"
    else:
        assert not "u" in opts

    #print "final options", opts

    bopts = opts.encode("ascii")
    copts = bopts
    # run the triangulation.
    with nogil:
        c_triangulate(copts, &tin, &tout, &tvorout)

    if "u" in opts:
        py_triunsuitable = NULL

    # extract all valid return arrays from the triangulateio structures,
    # wrap them in numpy arrays that free() the memory on their
    # deallocation (python garbage collection)
    # and fill the return dictionary with the arrays.
    ret = {}
    if not "N" in opts:
        ret["points"] = wrap(tout.pointlist, np.NPY_DOUBLE, 2,
                [tout.numberofpoints, 2])
    if not ("N" in opts or tout.numberofpointattributes == 0):
        ret["pointattributes"] = wrap(tout.pointattributelist,
                np.NPY_DOUBLE, 2, [tout.numberofpoints,
                    tout.numberofpointattributes])
    if not ("N" in opts or "B" in opts):
        ret["pointmarkers"] = wrap(tout.pointmarkerlist, np.NPY_INT, 1,
                [tout.numberofpoints])
    if not "E" in opts:
        ret["triangles"] = wrap(tout.trianglelist, np.NPY_INT, 2,
                [tout.numberoftriangles, tout.numberofcorners])
    if not "E" in opts and (
            tout.numberoftriangleattributes != 0 or "A" in opts):
        ret["triangleattributes"] = wrap(tout.triangleattributelist,
                np.NPY_DOUBLE, 2, [tout.numberoftriangles,
                    tout.numberoftriangleattributes])
    if "n" in opts:
        ret["neighbors"] = wrap(tout.neighborlist, np.NPY_INT, 2,
                [tout.numberoftriangles, 3])
    if ("p" in opts or "c" in opts) and not "P" in opts:
        ret["segments"] = wrap(tout.segmentlist, np.NPY_INT, 2,
                [tout.numberofsegments, 2])
    if ("p" in opts or "c" in opts) and not ("P" in opts or "B" in opts):
        ret["segmentmarkers"] = wrap(tout.segmentmarkerlist, np.NPY_INT,
                1, [tout.numberofsegments])
    if "e" in opts:
        ret["edges"] = wrap(tout.edgelist, np.NPY_INT, 2,
                [tout.numberofedges, 2])
    if "e" in opts and not "B" in opts:
        ret["edgemarkers"] = wrap(tout.edgemarkerlist, np.NPY_INT, 1,
                [tout.numberofedges])
    if "v" in opts:
        ret["voronipoints"] = wrap(tvorout.pointlist, np.NPY_DOUBLE, 2,
                [tvorout.numberofpoints, 2])
        if not ("N" in opts or tvorout.numberofpointattributes == 0):
            ret["voronipointattributes"] = wrap(tvorout.pointattributelist,
                    np.NPY_DOUBLE, 2, [tvorout.numberofpoints,
                        tvorout.numberofpointattributes])
        ret["voroniedges"] = wrap(tvorout.edgelist, np.NPY_INT, 2,
                [tvorout.numberofedges, 2])
        ret["voroninorms"] = wrap(tvorout.normlist, np.NPY_DOUBLE, 2,
                [tvorout.numberofedges, 2])

    return ret


# global....
@cython.boundscheck(False)
@cython.wraparound(False)
cdef public int triunsuitable(double* triorg, double* tridest, double* triapex,
        double area) with gil:
    """copied from triangle.c"""
    cdef double dxoa, dxda, dxod
    cdef double dyoa, dyda, dyod
    cdef double oalen, dalen, odlen
    cdef double maxlen

    if py_triunsuitable != NULL:
        return (<object>py_triunsuitable)(
                (triorg[0], triorg[1]), (tridest[0], tridest[1]), 
                (triapex[0], triapex[1]), area)

    dxoa = triorg[0] - triapex[0]
    dyoa = triorg[1] - triapex[1]
    dxda = tridest[0] - triapex[0]
    dyda = tridest[1] - triapex[1]
    dxod = triorg[0] - tridest[0]
    dyod = triorg[1] - tridest[1]
    # Find the squares of the lengths of the triangle's three edges.
    oalen = dxoa*dxoa + dyoa*dyoa
    dalen = dxda*dxda + dyda*dyda
    odlen = dxod*dxod + dyod*dyod
    # Find the square of the length of the longest edge.
    maxlen = max(dalen, oalen)
    maxlen = max(odlen, maxlen)

    if maxlen > 0.05*(triorg[0]*triorg[0] + triorg[1]*triorg[1])+0.02:
        return 1
    else:
        return 0
