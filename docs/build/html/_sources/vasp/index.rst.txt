==============
VASP RESOURCES
==============

.. toctree::
  :maxdepth: 2

  hardware_requirements

This section contains a variety of different notes on using the plane augmented wave (PAW) density functional theory (DFT) code called VASP.

Compilation of VASP

Running a VASP simulation

VASP cookbooks:
- energy cutoff convergence
- kpoint mesh convergence
- kpoint density per reciprocal atom
- structural minimization
- nudged elastic band calculation
- calculation of elastic constants
  
In general use, use IBRION=6, ISIF=3, POTIM=0.015in general use


David Holec, Martin Friák, Jörg Neugebauer, and Paul H. Mayrhofer. Phys. Rev. B 85, 064101
R., Zhu J., Ye H. Calculations of single-crystal elastic constants made simple. Comput. Phys. Commun. 2010;181(3):671–675
VASP checklist: This is a list of tasks which should be done to determine whether or not your VASP simulation worked.


A VASP simulation requires the following files at a minimum: INCAR, POSCAR, POTCAR, and KPOINTS.

