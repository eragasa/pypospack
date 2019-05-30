import yaml, copy, pathlib
import numpy as np
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.task.vasp as tsk_vasp
import os

if __name__ == '__main__':
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

    directory='MgO_NaCl'
    structure_filename=os.path.join(
            os.getcwd(),'rsrc','MgO_NaCl_prim.vasp')
    encut_conv = tsk_vasp.VaspEncutConvergence(
            directory=directory,
            structure=structure_filename,
            xc=xc,
            incar_dict=incar_dict,
            slurm_dict=slurm_dict,
            full_auto=full_auto)
