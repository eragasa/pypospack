import pytest


structure_filename = 'MgO_NaCl_prim.vasp'
xc = 'GGA'

full_auto=True
incar_dict = {}
incar_dict['ismear'] = 0
incar_dict['sigma'] = 0.05
incar_dict['ispin'] = 1

slurm_dict = {}
slurm_dict['email'] = 'eragasa@ufl.edu'
slurm_dict['qos'] = 'phillpot'
slurm_dict['ntasks'] = 16
slurm_dict['time'] = "1:00:00"


