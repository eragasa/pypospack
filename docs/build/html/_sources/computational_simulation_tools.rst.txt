##############################
Computational Simulation Tools
##############################

Contents
========

.. toctree::
   :maxdepth: 2

   crystallography/index
   dft/index
   eip/index
   md/index
   ld/index
   calc_material_properties/index
   pot_dev/index
   vasp/index
   lammps/index
   gulp/index

The purpose of these notes is to provide an introduction to the physical properties of solids, which are of extraordinarily important in the modern world.  It focuses upon fundamental, unifying concepts in solid-state physics but with a computational approach to under the properties of nuclei and electrons in solids rather than the typical experimental approach found in most solid state physics books.  These notes look to provide a practical introduction to computational materials to establish the basic principles, to describe phenomena responsible for the importance of solids in science and technology, while at the same time orienting the discussion toward understanding computational materials techniques and using them for explanatory purposes as well as to predict behavior.

The purpose of these notes is not intended to be exhaustive there exists many existing textbooks on computational materials as well the theory in materials science.  As a result, I have attempted to provide a fairly extensive bibliography which include recommended textbooks, as well as citation to the original journal articles.  Since this is an online reference, I intend to develop this website as a series of monographs and tutorials.

Introduction to Computational Tools
===================================

General Reference Books for Computational Material Science
----------------------------------------------------------
* :cite:`lee2016_book_comp_mse` Lee, J.G.  Computational Materials Science: an Introduction. 2016. CRC Press.

Crystallography
---------------

:ref:`The Crystal Lattice <crystal_lattice>`

:ref:`Reciprocal Lattice <reciprocal_lattice>`

For building crystal structres in a programatic way, the atomistic simulation environment is well-developed package.

Atomic Simulation Environment.  `Building Things <https://wiki.fysik.dtu.dk/ase/ase/build/build.html#common-bulk-crystals>`_


.. :ref:`Python excercises in crystallography <python_crystallography>`


Density Functional Theory
-------------------------
:ref:`Introduction to Density Functional Theory, by Eugene Ragasa <dft_intro>`

`VASP <https://www.vasp.at>`_

`pypospack VASP examples <pypospack/examples/vasp/index>`

Molecular Dynamics
------------------
:ref:`Introduction to Molecular Dynamics, by Eugene Ragasa <intro_md>`

`Introduction to Molecular Dynamics Simulation, by Micheal P. Allen <https://udel.edu/~arthij/MD.pdf>`_

`LAMMPS <http://lammps.sandia.gov>`_

Lattice Dynamics
----------------
The study of vibrations of atoms inside crystals, lattice dynamics, is basic to many fields of study in solid state physics and materials science.

`Introduction to Lattice Dynamics, by Martin T. Dove <https://www.amazon.com/Introduction-Lattice-Dynamics-Cambridge-Chemistry/dp/0521398940>`_ is a well written book accessible to advanced undergraduates, graduate students, and reserach workers looking to understand phonons from an approachable but fairly rigorous perspective.

`gulp <https://gulp.curtin.edu.au/gulp/>`_ is software written in Fortran which can do a variety of molecular dynamics type tasks, but really excells in calculated phonons and phonon density of states for a variety of empirical interatomic potentials.

`phonopy <https://atztogo.github.io/phonopy/>`_

Material Properties and Defects in Crystals
===========================================

These links are useful for building defect stuctures.

* Atomic Simulation Environment.  `General crystal structure <https://wiki.fysik.dtu.dk/ase/ase/lattice.html>`_

.. * :ref:`Structural minimization <e_min>`

.. * :ref:`Point defects <p_defects>`

.. * :ref:`Calculating Surface Energies <calc_surface_energy>`

.. * :ref:`Nudged Elastic Band <neb>`

Development of Interatomic Potentials
=====================================

.. * :ref:`Typical Approach to Potential Parameterization`

.. * :ref:`Pareto Approach to Potential Parameterization``

References
==========

.. bibliography:: computational_simulation_tools.bib
