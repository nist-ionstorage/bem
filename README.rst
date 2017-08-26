BEM
===


License
-------

Triangle (in triangle/)
.......................

These programs may be freely redistributed under the condition that the
copyright notices (including the copy of this notice in the code
comments and the copyright notice printed when the `-h` switch is
selected) are not removed, and no compensation is received.  Private,
research, and institutional use is free.  You may distribute modified
versions of this code UNDER THE CONDITION THAT THIS CODE AND ANY
MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE
ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE
WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS.
Distribution of this code as part of a commercial system is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR.  (If you are not directly
supplying this code to a customer, and you are instead telling them how
they can obtain it for free, then you are not required to make any
arrangement with me.)

Fastlap (in fastlap/)
.....................

This software is being provided to you, the LICENSEE, by the Massachusetts
Institute of Technology (M.I.T.) under the following license. By
obtaining, using and/or copying this software, you agree that you have
read, understood, and will comply with these terms and conditions:

Permission to use, copy, modify and distribute this software and its
documentation for any purpose and without fee or royalty is hereby granted,
provided that you agree to comply with the following copyright notice and
statements, including the disclaimer, and that the same appear on ALL
copies of the software and documentation, including modifications that you
make for internal use or for distribution:

Copyright 1992 by the Massachusetts Institute of Technology. All rights
reserved.

THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES, EXPRESS OR IMPLIED. By way of example, but not limitation,
M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS
FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR
DOCUMENTATION WILL NOT INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS,
TRADEMARKS OR OTHER RIGHTS.

The name of the Massachusetts Institute of Technology or M.I.T. may NOT
be used in advertising or publicity pertaining to distribution of the
software. Title to copyright in this software and any associated
documentation shall at all times remain with M.I.T., and USER agrees to
preserve same.

Written by: K. Nabors, T. Korsmeyer, and J. White

Code written by NIST employees (examples/, inventor/)
.....................................................

This software was developed at the National Institute of Standards and
Technology (NIST) by employees of the Federal Government in the course
of their official duties. Pursuant to title 17 Section 105 of the United
States Code, this software is not subject to copyright protection and is
in the public domain. NIST assumes no responsibility whatsoever for its
use by other parties, and makes no guarantees, expressed or implied,
about its quality, reliability, or any other characteristic.

Wrappers, integration and driver layers (bem/)
..............................................

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Introduction
------------

Wraps the fastlap BEM (boundary element method aka BIE boundary integral
equations) + FMM (fast multipole method) + Laplace kernel implementation
from http://www.rle.mit.edu/cpg/research_codes.htm The original code is
in doc/fastlap_fl-2.0-22oct96 for reference.


General Notes
-------------

* The revised workflow for 3D electrostatics is currently in the
  "big_rewrite" branch in svn. It may move to trunk later.

* Most of the bugs and problems mentioned in "BEM software -
  Development" on the wiki are solved.

* The old .cpx "init" style configuration file is gone. Replaced with
  Python code. More flexible, less coding overhead when implementing new
  features. Compare the old cpx with the new python example.

* The input geometry is not the inventor-meshed .cpy file anymore but
  the face loops .ele or .trap or inventor exported .stl. You need to
  re-export your files or write a cpy-to-trap convertor (should not be too
  difficult).

* The meshing is now done preferably with Jonathan Shewchuk's triangle
  code. It is integrated into the Python package and compiled into a
  python module, see triangle and pytriangle.pyx

* The Fastlap code is also integrated and compiled into a python module
  (see fastlap/ and fastlap.pyx and called directly without exporting and
  importing intermediate file formats.

* The meshing can be adaptive: for a given potential configuration, the
  charge distribution for a small initial mesh is solved and the mesh is
  refined such that each triangle in the final mesh contributes equally to
  the field at the observation point (since the triangles are also
  constrained by the geometry, min angle etc, the final number of
  triangles is typically larger than the desired one).

* Triangle areas can also be constrained via differently shaped constraint
  fields (Box, Sphere, Cylinder or BYO).


STL
---

Ryan Bowler, 2014

If you wish to generate your own STL for the SimpleTrap, there are some
important features. The code scales the trap as if it is in units if
microns, so when exporting to STL in Inventor, be sure to choose microns
for the Units. There are no curved surfaces, so Surface and Normal
deviations are not important. Set max edge length to something
reasonable or else triangles are too small or large (40 microns, the ion
height from the surface, is a good choice). Choose a low aspect ratio
and make sure to export the colors.


Notes
-----

Old notes regarding multipole expansion and jumps in potentials/fields:

    The 'slfcc - Precise.exe' version is meant to solve the following
    problem. It can happen that the center of your trapping region is
    right on the boundary between two "cells" of the tree structure
    built by FastLap for the multipole-accelerated algorithm. In this
    case the calculated potentials and fields will show tiny "jumps" in
    their values when going across this boundary. This has usually no
    noticeable effect on potentials, but can be noticeable on the field
    and hence pseudopotential.

    One way to solve this problem is to add dummy electrodes on the side
    of your real electrodes, so that the spatial structure of the tree
    is shifted a bit. This would displace the cell boundary out of the
    center of your trap.

    The other way is using 'slfcc - Precise.exe', which skips the
    multipole acceleration procedure when calculating potentials and
    fields. In other words, it does an exact calculation based on the
    solved charge distribution, without using any tree structure. This
    increases the computing time and memory requirements, but yields a
    slightly more precise result. Note that the charge solving part of
    the algorithm is not modified (= it uses multipole acceleration,
    with a depth set in script 'runBEM.py').

    -> In the new python code this is achieved by passing "num_lev=1" to
    Job.simulate().


Some File descriptions
----------------------

The cpx&cpy.reg file assumes a root directory C:\BEMcode
The vtk.reg file assumes a directory C:\Program Files\ParaView\

Examples\TesSphere_1mm\       1mm tessellated sphere
Examples\SimpleTrap\          Simple Signe style trap
Examples\Skull trap\          Skull trap outline to test Inventor import macros
