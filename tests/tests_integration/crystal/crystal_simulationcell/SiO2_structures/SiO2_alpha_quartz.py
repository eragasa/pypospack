
# alpha-quartz
# lattice: trigonal
# sg: P3_2 21
# sg_num: 154
# a=4.9134
# c=5.4052
# alpha = 90
# beta = 90
# gamma = 120

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp

filename = 'xuexiong_structures/SiO2_a_quartz.vasp'

poscar = vasp.Poscar()
poscar.read(filename)


# atomic basis
# Si 3a 0.4699 0.0000 0.0000 
# O1 6c 0.4141 0.2681 0.1188 
# beta-quartz
# lattice: 
