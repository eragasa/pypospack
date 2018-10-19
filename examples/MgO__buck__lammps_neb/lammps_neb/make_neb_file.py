#!/bin/env python
import sys

# this file will only work for this NEB simulation and has not been written
# generally enough to work for any NEB simulation
# eragasa @ 2017/02/03

# read line
fname = str(sys.argv[1])
f = open(fname,'r')
lines = f.readlines()
f.close()

lines = [l.strip() for l in lines]

# get number of atoms 
n_atoms = lines[2].split(' ')[0]
n_atoms = int(n_atoms)
begin_line = 17

# get atomic position
p_atoms = []
for i in range(begin_line,begin_line+n_atoms):
    line = lines[i].split(' ')
    atom_id = int(line[0])    
    x = float(line[3])
    y = float(line[4])
    z = float(line[5])
    p_atoms.append([atom_id, x, y, z])


str_out = str(n_atoms) + "\n"
for p in p_atoms:
    str_out += " ".join([str(v) for v in p]) + "\n"

f = open(fname+".neb",'w')
f.write(str_out)
f.close()
