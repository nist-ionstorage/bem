#!/usr/bin/python
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


import logging, numpy as np, argparse
from multiprocessing import Pool

from bem import Grid, Simulation, Result


def evaluate(prefix, electrodes):
    logging.basicConfig(level=logging.DEBUG,
        format="%(asctime)s %(name)s %(processName)s[%(process)s] "
               "%(filename)s:%(funcName)s() %(levelname)s: %(message)s",
        filename="%s.log" % prefix, filemode="w")
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter("%(processName)s[%(process)s] "
        "%(filename)s:%(funcName)s(): %(message)s"))
    logging.getLogger().addHandler(sh)

    n, s = 2*10, .1 # FIXME: parametrize
    grid = Grid(center=(0, 0, 1.5), step=(s, s, s), shape=(n, n, n))

    pool = Pool()
    c = map(run_job,
            ((grid, prefix, (e, True, True, True)) for e in electrodes))

    for e in electrodes:
        Result.view(prefix, e+"_new")

def run_job(args):
    grid, prefix, job = args
    logging.info("starting job %s", job)

    result = Result.from_vtk(prefix, job[0])
    sim = Simulation(result.mesh)
    result = sim.simulate(job, grid, sing=result.mesh_charge,
            fmm_on_eval=False)
    result.name += "_new" # FIXME: parametrize
    result.to_vtk(prefix)

def main():
    p = argparse.ArgumentParser(
            description="evaluate a bem simulation result (charge"
            " distribution) in another grid")
    p.add_argument("prefix", help="data path prefix")
    p.add_argument("electrodes", nargs="+", help="electrode name")
    args = p.parse_args()
    evaluate(args.prefix, args.electrodes)

if __name__ == "__main__":
    main()
    from pyface.api import GUI
    GUI().start_event_loop()
