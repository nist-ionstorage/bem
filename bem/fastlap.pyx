# -*- coding: utf8 -*-
#
#   bem: triangulation and fmm/bem electrostatics tools 
#   pyfastlap - python bindings for fastlap
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

"""BEM FMM Laplace solver"""

import sys
import cython

import numpy as np
cimport numpy as np

np.import_array()


cdef extern from "fastlap.h":
    cdef int c_fastlap "fastlap"(int *plhsSize, int *prhsSize,
            int *pnumSing, double *px, int *pshape, int *pdtype,
            int *plhsType, int *prhsType,
            int *plhsIndex, int *prhsIndex,
            double *plhsVect, double *prhsVect,
            double *pxf, double *pxnrm,
            int *pnumLev, int *pnumMom, int *pmaxItr,
            double *ptol, int *pjob, double* pAreas) nogil

    cdef enum PANEL_TYPES:
        FL_POINT "POINT" = 1
        FL_TRIANGLE "TRIANGLE" = 3
        FL_QUADRILAT "QUADRILAT" = 4

    cdef enum SOURCE_TYPES:
        FL_NO_SOURCE "NO_SOURCE" = 0
        FL_POINT_SOURCE "POINT_SOURCE" = 1
        FL_CONSTANT_SOURCE "CONSTANT_SOURCE" = 11
        FL_CONSTANT_DIPOLE "CONSTANT_DIPOLE" = 12
        FL_LINEAR_SOURCE "LINEAR_SOURCE" = 21
        FL_LINEAR_DIPOLE "LINEAR_DIPOLE" = 22

    cdef enum JOB_TYPES:
        FL_FIELD "FIELD" = 0
        FL_GREEN "GREEN" = 1
        FL_INDIRECT "INDIRECT" = 2

cdef extern from "memtracker.h":
    cdef void mtinit() nogil
    cdef void mtclear() nogil

cdef extern from "fastlap_support.h":
    cdef void Dcentroid(int shape, double *pc, double *xcout) nogil


POINT            = FL_POINT
TRIANGLE         = FL_TRIANGLE
QUADRILAT        = FL_QUADRILAT
                
NO_SOURCE        = FL_NO_SOURCE
POINT_SOURCE     = FL_POINT_SOURCE
CONSTANT_SOURCE  = FL_CONSTANT_SOURCE
CONSTANT_DIPOLE  = FL_CONSTANT_DIPOLE
LINEAR_SOURCE    = FL_LINEAR_SOURCE
LINEAR_DIPOLE    = FL_LINEAR_DIPOLE
                
FIELD            = FL_FIELD
GREEN            = FL_GREEN
INDIRECT         = FL_INDIRECT


ctypedef int intc_t # np.int_t is a long, int (C) is 32 bits

@cython.boundscheck(False)
@cython.wraparound(False)
def fastlap(
    np.ndarray[np.double_t, ndim=3, mode="c"] x not None, # M, 4, 3
    np.ndarray[intc_t, ndim=1, mode="c"] shape not None, # M
    np.ndarray[intc_t, ndim=1, mode="c"] lhs_type not None, # M
    np.ndarray[intc_t, ndim=1, mode="c"] rhs_type not None, # M
    np.ndarray[np.double_t, ndim=1, mode="c"] rhs_vect not None, # <=M*n, n=1
    np.ndarray[np.double_t, ndim=2, mode="c"] xf not None, # N, 3
    np.ndarray[intc_t, ndim=1, mode="c"] type not None, # N
    int job,
    np.ndarray[np.double_t, ndim=2, mode="c"] xnrm = None, # N, 3
    np.ndarray[intc_t, ndim=1, mode="c"] lhs_index = None, # M*n, n=1
    np.ndarray[intc_t, ndim=1, mode="c"] rhs_index = None, # M*n, n=1
    np.ndarray[np.double_t, ndim=1, mode="c"] lhs_vect = None, # N
    int ret_areas = False,
    int num_mom = 2, int num_lev = 4,
    int max_iter = 32, double tol = 1e-4):
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
    cdef int num_iter
    cdef int M = x.shape[0], N = xf.shape[0]
    cdef int lhs_size = N, rhs_size = rhs_vect.shape[0], num_sing = M
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] areas = None
    cdef double *xnrm_ = NULL, *areas_ = NULL

    if lhs_index is None:
        lhs_index = np.arange(M, dtype=np.intc)
    if rhs_index is None:
        rhs_index = np.arange(M, dtype=np.intc)
    if lhs_vect is None:
        lhs_vect = np.empty((N,), dtype=np.double)
    if xnrm is not None:
        xnrm_ = &xnrm[0, 0]
    if ret_areas:
        areas = np.empty((M,), dtype=np.double)
        areas_ = &areas[0]

    assert x.shape[0] == M
    assert x.shape[1] == 4
    assert x.shape[2] == 3
    assert shape.shape[0] == M
    assert lhs_type.shape[0] == M
    assert rhs_type.shape[0] == M
    assert lhs_index.shape[0] == M
    assert rhs_index.shape[0] == M
    assert rhs_vect.shape[0] <= M
    assert lhs_vect.shape[0] == N
    assert xf.shape[0] == N
    assert xf.shape[1] == 3
    assert type.shape[0] == N

    with nogil:
        mtinit()
        num_iter = c_fastlap(&lhs_size, &rhs_size, &num_sing,
                &x[0, 0, 0], &shape[0], &type[0],
                &lhs_type[0], &rhs_type[0],
                &lhs_index[0], &rhs_index[0],
                &lhs_vect[0], &rhs_vect[0], &xf[0, 0], xnrm_,
                &num_lev, &num_mom, &max_iter, &tol, &job, areas_)
        mtclear()
    return num_iter, tol, lhs_vect, areas


@cython.boundscheck(False)
@cython.wraparound(False)
def centroid(
        np.ndarray[np.double_t, ndim=3, mode="c"] pc not None,
        np.ndarray[intc_t, ndim=1, mode="c"] shape not None,
        np.ndarray[np.double_t, ndim=2, mode="c"] xc = None
        ):
    assert shape.shape[0] == pc.shape[0]
    if xc is None:
        xc = np.empty((pc.shape[0], 3), dtype=np.double)
    assert xc.shape[0] == pc.shape[0]
    with nogil:
        for i in range(pc.shape[0]):
            Dcentroid(shape[i], &pc[i, 0, 0], &xc[i, 0])
    return xc
