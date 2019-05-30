====================
Tutorials for mpi4py
====================

The tutorials written here are adapted from MPI tutorials written in different different programming languages, but then adapted here for provide an mpi4py tutorial.  The purpose of writing these pages was for me to learn both MPI and mpi4py.  So these pages are notes on what I learned as I learned it.

Installation on OSX
===================

Installing mpi4py on anaconda is fairly simple.  First ensure that you have don't have any MPI stuff installed.

conda list | grep mpi

conda install mpi4py

Installation on a Linux Cluster
===============================

The scripts here include SGE and slurm scripts.


.. toctree::
   :maxdepth: 1

   mpi4py_01
