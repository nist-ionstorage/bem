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


import argparse
from pyface.api import GUI
from bem import Result


def main():
    p = argparse.ArgumentParser(description="visualize bem simulation result")
    p.add_argument("prefix", help="data path prefix")
    p.add_argument("electrodes", nargs="+", help="electrode name")
    args = p.parse_args()
    for e in args.electrodes:
        Result.view(args.prefix, e)
    GUI().start_event_loop()


if __name__ == "__main__":
    main()
