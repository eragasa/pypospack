=====================
An Example using VASP
=====================

The Vienna Ab Initio Simulation Package (VASP) is a computer program for electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.  VASP computes an approximate solution to the many-body Schrodinger's equation, either with density functional theory (DFT) by solving the Kohn-Sham equations, or within the Hartree-Fock (HF) approximation by solving the Roothaan equations.  Writing a software package which completely encapsulates VASP is likely to be brittle and will likely implement many features which are not necessary.

This software development process first starts, by implementing the simplest calculation the software will support with the minimum standard tags.  A minimal VASP simulation requires a small number of files: POSCAR, POTCAR, INCAR, and KPOINTS file.

Running the Initial Simulation
==============================

Let us look at a typical POSCAR file for VASP by looking at a magnesium oxide (MgO), the ground state of MgO.  The first line is a comment line which we could possibly use to encode information about the structure, which typically cannot be simply determined by inspection of the file.  We separate the POSCAR file into the following components.
