..
    .. image:: https://travis-ci.org/farrelljd/pele.svg?branch=modernise
        :target: https://travis-ci.org/pele-python/pele?branch=master

    .. image:: https://coveralls.io/repos/pele-python/pele/badge.png?branch=master
        :target: https://coveralls.io/r/pele-python/pele?branch=master

    .. image:: https://landscape.io/github/pele-python/pele/master/landscape.svg
       :target: https://landscape.io/github/pele-python/pele/master
       :alt: Code Health

pele : Python Energy Landscape Explorer
+++++++++++++++++++++++++++++++++++++++

Tools for global optimization and energy landscape exploration.

Source code: https://github.com/pele-python/pele

Documentation: http://pele-python.github.io/pele/

.. figure:: ./images/lj38_gmin_dgraph.png

  Images: The global minimum energy structure of a 38 atom Lennard-Jones cluster.  On
  the right is a disconnectivity graph showing a visualization of the energy
  landscape.  The competing low energy basins are shown in color.

pele is a python partial-rewriting of GMIN, OPTIM, and PATHSAMPLE: fortran
programs written by David Wales of Cambridge University and collaborators
(http://www-wales.ch.cam.ac.uk/software.html).  

Description
===========
pele has tools for energy minimization, global optimization, saddle point
(transition state) search, data analysis, visualization and much more.  Some of
the algorithms implemented are:

1. Basinhopping global optimization

#. Energy minimization

   - LBFGS
   - FIRE

#. Single ended saddle point search

   - Hybrid Eigenvector Following
   - Dimer method

#. Double ended saddle point search

   - Nudged Elastic Band (NEB)
   - Doubly Nudged Elastic Band (DNEB)

#. Disconnectivity Graph visualization

#. Structure alignment algorithms

#. Thermodynamics (e.g. heat capacity) via the Harmonic Superposition Approximation

#. Transition rates analysis

Installation
============

Building and installing pele requires the meson build system, C++ and FORTRAN compilers for building compiled
extensions, and a Python installation (version >=3.9, <3.14).
These dependencies can be satisfied by creating a new Conda environment from the packaged ``environment.yml`` file::

    conda env create -f environment.yml
    conda activate pele

With the dependencies installed, installing pele is as simple as::

    git clone https://github.com/pele-python/pele.git
    cd pele
    pip install .

Options
--------

To make tests available, replace ``.`` with ``.[testing]``, *i.e.*::

    pip install .[testing]

and then, outside of this source tree, run::

    nosetests pele

We also have tests for our c++ code writen in pure c++.  These are stored in
the folder ``cpp_tests/`` and can be compiled using CMake. These tests have not been tested recently.

To use the new LAMMPS interface::

    pip install .[lammps]

making the pure Python LAMMPS potential class and related system classes available.

To use the compiled LAMMPS interface::

    pip install .[lammps] -Csetup-args="-Dlammps=enabled"

making the compiled c++ LAMMPS potential class available.

Pip will download and install the LAMMPS Python package if one is not already on the path, along with MPICH headers.
If the LAMMPS python package has been compiled and installed from source, the relevant MPI headers must also be on the
path.

Running
=======

You can find examples of how to run pele in the examples folder.  More information can be found in the documentation at

http://pele-python.github.com/pele/


Notes
=====
In the long-long-ago, pele was known as pygmin.