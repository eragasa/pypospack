import pytest
import shutil, copy, os
import numpy as np
import pypospack.task.vasp as tsk_vasp
import pypospack.io.vasp as vasp
import pypospack.dft as dft
from collections import OrderedDict

def cleanup(directory):
    shutil.rmtree(directory)

# this approach is taken from the VASP methodology of generating automated kpoint meshes
# Ref: https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html

#xc = 'GGA'

#full_auto=False
#3incar_dict = {}
#incar_dict['ismear'] = 0
#incar_dict['sigma'] = 0.05
#incar_dict['ispin'] = 1

#slurm_dict = {}
#slurm_dict['email'] = 'eragasa@ufl.edu'
#slurm_dict['qos'] = 'phillpot'
#slurm_dict['ntasks'] = 16
#slurm_dict['time'] = "1:00:00"

# not sure if this works
#poscar = vasp.Poscar(structure_filename)


if __name__ == "__main__":
    directory = 'SiO2_aQ'
    structure_filename = os.path.join(
            'SiO2_structures',
            'SiO2_aq_cubic.vasp')

    poscar = vasp.Poscar()
    poscar.read(structure_filename)
    print('structure_filename:structure_filename')
    print('--- REAL SPACE LATTICE ---')
    print('\ta1:',poscar.H[0,:],'->',poscar.a1)
    print('\ta2:',poscar.H[1,:],'->',poscar.a2)
    print('\ta3:',poscar.H[2,:],'->',poscar.a3)
    print('--- RECIPROCAL SPACE LATTICE ---')
    print('\tb1:',poscar.b1,'->',np.dot(poscar.b1,poscar.b1)**0.5)
    print('\tb2:',poscar.b2,'->',np.dot(poscar.b2,poscar.b2)**0.5)
    print('\tb2:',poscar.b3,'->',np.dot(poscar.b3,poscar.b3)**0.5)
    print('--- KPOINT MESHES ---')

    # linear kpoint density
    kwargs = {
            'simulation_cell':poscar,
            'rho_min':1,
            'rho_max':10,
            'd_rho':0.1,
            'kpoint_min':3,
            'kpoint_max':15}
    kpoint_meshes = dft.determine_kpoint_meshes(**kwargs)
    # print kpoint_meshes
    for k,v in kpoint_meshes.items():
        print(k,':',v)
