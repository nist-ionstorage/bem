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

import numpy as np, logging
from .fastlap import fastlap, centroid, CONSTANT_SOURCE, NO_SOURCE, INDIRECT, FIELD, TRIANGLE


def sing_pot_field(x, pot=None, xm=None, shape=None, sing=None, segsize=None,
        do_field=True, fmm_on_eval=True, num_mom=4, num_lev=4,
        max_iter=500, tol=1e-5):
    """
    Calculates the charges on the panels at `x` with shape `shape` given
    their potentials `pot`. Then calculates the potential and field at
    points `xm`.

    Returns singularity strenghts (charge densities) and areas of the panels,
    then potentials and field components at `xm`
    """

    x = np.ascontiguousarray(x)
    m = x.shape[0] # 20e3 takes 4GB, 3min

    if shape is None:
        shape = TRIANGLE*np.ones((m,), dtype=np.intc)
    xf = centroid(x, shape)
    constant_source = CONSTANT_SOURCE*np.ones((m,), dtype=np.intc)
    no_source = NO_SOURCE*np.ones((m,), dtype=np.intc)

    index = np.arange(m, dtype=np.intc)

    opts = dict(x=x, shape=shape, lhs_index=index, rhs_index=index,
            num_mom=num_mom, num_lev=num_lev, max_iter=max_iter, tol=tol)

    if pot is not None:
        pot = np.ascontiguousarray(pot)

        n_itr, tol, sing, areas = fastlap(
                lhs_type=constant_source, rhs_type=no_source,
                rhs_vect=pot, xf=xf,
                type=np.zeros((m,), dtype=np.intc),
                job=INDIRECT, ret_areas=True, **opts)
        logging.info("%i triangles, gmres %i iters, %g tolerance", m, n_itr, tol)

        if xm is None:
            return sing, areas, xf
    else:
        sing = np.ascontiguousarray(sing)
        areas = None

    xm = np.ascontiguousarray(xm)
    n = xm.shape[0] # m*n/(20e3*25e3) takes 8GB, 3min

    if segsize is None:
        segsize = int(25e7/m) # about 4GB
    segsize = min(n, segsize)

    pots = np.empty((n,), dtype=np.double)
    typep = np.zeros((segsize,), dtype=np.intc)

    if do_field:
        field = np.empty((3, n), dtype=np.double)
        typef = np.ones((segsize,), dtype=np.intc)
        xnrm = np.zeros((segsize, 3), dtype=np.double)
        xnrm[:, 0] = 1.
    else:
        field = None

    if not fmm_on_eval:
        opts["num_lev"] = 1
    
    for i in range(0, n, segsize):
        segsize = min(segsize, n-i)
        n_itr, tol, _, _ = fastlap(
                lhs_type=no_source, rhs_type=constant_source,
                rhs_vect=sing, xf=xm[i:i+segsize],
                lhs_vect=pots[i:i+segsize],
                type=typep[:segsize], job=FIELD, **opts)
        logging.info("fastlap potential points: %i-%i/%i",
                i, i+segsize, n)
        if not do_field:
            continue
        for j in range(3):
            n_itr, tol, _, _ = fastlap(
                    lhs_type=no_source, rhs_type=constant_source,
                    rhs_vect=sing, xf=xm[i:i+segsize],
                    lhs_vect=field[j, i:i+segsize],
                    xnrm=np.roll(xnrm, j, axis=1)[:segsize],
                    type=typef[:segsize], job=FIELD, **opts)
            logging.info("fastlap field[%i] points: %i-%i/%i",
                    j, i, i+segsize, n)
    return sing, areas, xf, pots, field



