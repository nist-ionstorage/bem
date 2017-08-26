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


import logging, re, operator, os

import numpy as np
try:
    from tvtk.api import tvtk
    from mayavi.modules.surface import Surface
    from mayavi.modules.iso_surface import IsoSurface
except ImportError:
    import warnings
    warnings.warn("not tvtk found")

from .fastlap import (fastlap, centroid,
        TRIANGLE, NO_SOURCE, CONSTANT_SOURCE, INDIRECT, FIELD)

from .triangulation import Mesh
from .grid import Grid


class Result(object):
    configuration = None
    grid = None
    potential = None
    field = None
    pseudo_potential = None

    def to_vtk(self, prefix):
        """
        export the result to vtk
        the mesh and its triangle data go to prefix_name_mesh.vtk,
        the potential/field/pseudopotential go to prefix_name.vtk
        all arrays are named
        """
        if self.configuration is not None:
            self.configuration.to_vtk(prefix)
        sp = tvtk.StructuredPoints(
                origin=self.grid.get_origin(),
                spacing=self.grid.step,
                dimensions=self.grid.shape)
        spw = tvtk.StructuredPointsWriter(input=sp)
        spw.file_name = "%s_%s.vtk" % (prefix, self.configuration.name)
        #xidr = tvtk.XMLImageDataWriter(input=sp)
        for data_name in "potential field pseudo_potential".split():
            data = getattr(self, data_name)
            if data is None:
                continue
            if len(data.shape) > 3:
                data = data.T.reshape(-1, data.shape[0])
            else:
                data = data.T.flatten()
            d = tvtk.DoubleArray(name=data_name)
            d.from_array(data)
            sp.point_data.add_array(d)
        spw.write()
        logging.info("wrote %s", spw.file_name)

    @classmethod
    def from_vtk(cls, prefix, name):
        """
        read in a result object from available vtk files
        """
        obj = cls()
        obj.configuration = Configuration.from_vtk(prefix, name)
        spr = tvtk.StructuredPointsReader(
            file_name="%s_%s.vtk" % (prefix, name))
        spr.update()
        sp = spr.output
        for i in range(sp.point_data.number_of_arrays):
            da = sp.point_data.get_array(i)
            name = da.name
            data = da.to_array()
            step = sp.spacing
            origin = sp.origin
            shape = tuple(sp.dimensions)
            dim = shape[::-1]
            if da.number_of_components > 1:
                dim += (da.number_of_components,)
            data = data.reshape(dim).T
            setattr(obj, name, data)
        # FIXME: only uses last array's data for grid
        obj.grid = Grid(step, shape, origin+(np.array(shape)-1)/2.*step)
        return obj

    @staticmethod
    def view(prefix, name):
        """
        construct a generic visualization of base mesh, refined mesh,
        and potential/field/pseudopotential data in mayavi2
        requires running within mayavi2 or ipython with threads or
        something like::

            from pyface.api import GUI
            GUI().start_event_loop()

        in your script to interact with it.

        this is from a simplified macro recorded in mayavi2
        """
        try:
            engine = mayavi.engine
        except NameError:
            from mayavi.api import Engine
            engine = Engine()
            engine.start()

        if len(engine.scenes) == 0:
            engine.new_scene()

        scene = engine.scenes[0]

        base_mesh_name = "%s_mesh.vtk" % prefix
        if os.access(base_mesh_name, os.R_OK):
            base_mesh = engine.open(base_mesh_name)
            surface = Surface()
            engine.add_filter(surface, base_mesh)
            surface.actor.property.representation = 'wireframe'
            surface.actor.property.line_width = 1

        mesh_name = "%s_%s_mesh.vtk" % (prefix, name)
        if os.access(mesh_name, os.R_OK):
            mesh = engine.open(mesh_name)
            mesh.cell_scalars_name = 'charge'
            surface = Surface()
            engine.add_filter(surface, mesh)
            module_manager = mesh.children[0]
            module_manager.scalar_lut_manager.lut_mode = 'RdBu'
            module_manager.scalar_lut_manager.use_default_range = False
            r = np.fabs(module_manager.scalar_lut_manager.data_range).max()
            module_manager.scalar_lut_manager.data_range = [-r, r]
            surface.actor.property.backface_culling = True

        data_name = "%s_%s.vtk" % (prefix, name)
        if os.access(data_name, os.R_OK):
            data = engine.open(data_name)
            if "pseudo_potential" in data._point_scalars_list:
                data.point_scalars_name = "pseudo_potential"
            else:
                data.point_scalars_name = "potential"
            iso_surface = IsoSurface()
            engine.add_filter(iso_surface, data)
            module_manager = data.children[0]
            module_manager.scalar_lut_manager.lut_mode = 'Greys'
            iso_surface.contour.auto_contours = True
            iso_surface.contour.number_of_contours = 5
            try:
                iso_surface.contour.maximum_contour = 1e-2
            except:
                pass

        scene.scene.isometric_view()
        scene.scene.render()


class Configuration(object):
    """
    a simulation configuration: simulate given mesh for certain
    potentials on certain lecetrodes
    """
    mesh = None
    potentials = None
    name = None

    charge = None

    opts = None
    data = None

    def __init__(self, mesh, potentials, name=None):
        self.mesh = mesh
        self.potentials = potentials
        self.name = name

    @classmethod
    def select(cls, mesh, *electrodes):
        """
        yields unit-potential simulation configurations given regexps for
        electrode names.
        """
        names = mesh.keys() # needs to be an OrderedDict
        for i, name in enumerate(names):
            match = any(re.match(p, name) for p in electrodes)
            if match or not electrodes:
                potentials = np.zeros(len(names))
                potentials[i] = 1.
                obj = cls(mesh, potentials, name)
                yield obj

    def set_data(self):
        """
        prepare arrays to be passed to fastlap
        """
        x = np.ascontiguousarray(self.mesh.fastlap_points())
        m = x.shape[0]
        panel_potential = self.potentials[self.mesh.groups]
        shape = TRIANGLE*np.ones((m,), dtype=np.intc)
        panel_centroid = centroid(x, shape)
        constant_source = CONSTANT_SOURCE*np.ones((m,), dtype=np.intc)
        no_source = NO_SOURCE*np.ones((m,), dtype=np.intc)
        index = np.arange(m, dtype=np.intc)
        typep = np.zeros((m,), dtype=np.intc)
        self.opts = dict(x=x, shape=shape, lhs_index=index, rhs_index=index)
        self.data = dict(constant_source=constant_source,
                no_source=no_source, typep=typep,
                potential=panel_potential,
                centroid=panel_centroid)

    def solve_singularities(self, num_mom=4, num_lev=3, max_iter=200,
            tol=1e-5, **fastlap_opts):
        """
        solve singularties' strengths (charge) on the panels. For other
        options, see fastlap().
        """
        logging.info("solving singularities %s", self.name or
                self.potentials)
        self.set_data()
        fastlap_opts.update(self.opts)
        fastlap_opts.update(dict(num_mom=num_mom, num_lev=num_lev,
            max_iter=max_iter, tol=tol))
        n_itr, tol, self.charge, self.data["area"] = fastlap(
                lhs_type=self.data["constant_source"],
                rhs_type=self.data["no_source"],
                rhs_vect=self.data["potential"],
                xf=self.data["centroid"],
                type=self.data["typep"],
                job=INDIRECT, ret_areas=True, **fastlap_opts)
        logging.info("n_tri=%i, iter=%i, tol=%g",
                len(self.charge), n_itr, tol)
        return n_itr, tol

    def collect_charges(self):
        """
        returns total accumulated charge per electrode
        Q_i=C_ij*U_j
        if U_j = 1V*delta_ij, this is the column Ci of the capacitance matrix
        in CGS units, multiply by 4*pi*epsilon_0*length_scale to get SI Farad
        (or length_scale*1e7/c**2)
        """
        return np.bincount(self.mesh.groups,
                self.charge * self.data["area"])

    def adapt_mesh(self, triangles=1e3, opts="qQ", min_area=1e-4,
            max_area=1e4, observation=np.array([0., 0., 1.])):
        """
        refine the mesh such that triangles are sized according
        to the field strength they create at `observation`.
        About `triangles` triangles are created (typically more).
        Refinement is globally limited to between min_area and
        max_area per triangle.
        """
        self.solve_singularities(num_mom=2, tol=1e-2)
        distance2 = np.square(self.data["centroid"]-observation).sum(axis=-1)
        weight = np.absolute(self.charge)/distance2
        max_areas = (self.data["area"]*weight).sum()/(triangles*weight)
        np.clip(max_areas, min_area, max_area, out=max_areas)
        logging.info("estimate at least %i triangles",
                np.ceil(self.data["area"]/max_areas).sum())
        self.mesh.set_max_areas(max_areas)
        self.mesh = self.mesh.triangulate(opts=opts, new=True)

    def evaluate(self, xm, segsize=None, derivative=False,
            num_mom=4, num_lev=3, **fastlap_opts):
        """
        evaluate potential (derivative=False) or field (derivative=True)
        at the given points `xm`. Fragment into `segsize` points. For other
        options, see fastlap().
        """
        fastlap_opts.update(self.opts)
        fastlap_opts.update(dict(num_mom=num_mom, num_lev=num_lev))
        xm = np.ascontiguousarray(xm)
        n = xm.shape[0] # m*n/(20e3*25e3) takes 8GB, 3min

        if segsize is None:
            segsize = int(25e7/len(self.charge)) # about 4GB
        segsize = min(n, segsize)

        if derivative:
            field = np.empty((3, n), dtype=np.double)
            type = np.ones((segsize,), dtype=np.intc)
            xnrm_ = np.zeros((segsize, 3), dtype=np.double)
            xnrm_[:, 0] = 1
        else:
            field = np.empty((1, n), dtype=np.double)
            type = np.zeros((segsize,), dtype=np.intc)
            xnrm = None

        for i in range(0, n, segsize):
            segsize = min(segsize, n-i)
            for j in range(field.shape[0]):
                if derivative:
                    xnrm = np.roll(xnrm_, j, axis=1)[:segsize]
                fastlap(lhs_type=self.data["no_source"],
                        rhs_type=self.data["constant_source"],
                        rhs_vect=self.charge,
                        xf=xm[i:i+segsize],
                        lhs_vect=field[j, i:i+segsize],
                        xnrm=xnrm,
                        type=type[:segsize],
                        job=FIELD, **fastlap_opts)
                logging.info("derivative=%s, j=%i: %i-%i/%i",
                        derivative, j, i, i+segsize, n)
        if not derivative:
            field = field[0]
        return field

    def simulate(self, grid, potential=True, field=False,
            pseudopotential=True, **fastlap_opts):
        """performs the simulations for observation
        points in `grid`. other options passed
        down to `fastlap()`.
        returns a `Result()` object containing the data.
        use num_lev=1 for fmm_on_eval=False
        """

        res = Result()
        res.configuration = self
        res.grid = grid

        xm = grid.to_points()

        if potential:
            pot = self.evaluate(xm, derivative=False, **fastlap_opts)
            res.potential = pot.reshape(grid.shape)
        if field:
            field = self.evaluate(xm, derivative=True, **fastlap_opts)
            res.field = field.reshape((3,) + grid.shape)
            if pseudopotential:
                pp = np.square(field).sum(axis=0).reshape(grid.shape)
                res.pseudo_potential = pp

        logging.info("done with job %s, %s, %s, %s", self.name, potential, field,
                pseudopotential)
        return res

    def to_vtk(self, prefix):
        self.mesh.to_vtk("%s_%s" % (prefix, self.name),
            area=self.data["area"],
            potential=self.data["potential"],
            charge=self.charge)

    @classmethod
    def from_vtk(cls, prefix, name):
        mesh, datasets = Mesh.from_vtk("%s_%s" % (prefix, name))
        pot = datasets.get("potential")
        if pot is not None:
            eles, reles = np.unique(mesh.groups, return_index=True)
            potentials = pot[reles]
        else:
            potentials = None
        obj = cls(mesh, potentials, name)
        obj.data = dict(area=datasets.get("area"))
        obj.charge = datasets.get("charge")
