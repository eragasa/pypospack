import pytest
import shutil, copy
import numpy as np
import pypospack.task.vasp as tsk_vasp
import pypospack.io.vasp as vasp

def cleanup(directory):
    shutil.rmtree(directory)

directory = 'MgO_NaCl'
structure_filename = '../rsrc/MgO_NaCl_prim.vasp'
xc = 'GGA'

full_auto=False
incar_dict = {}
incar_dict['ismear'] = 0
incar_dict['sigma'] = 0.05
incar_dict['ispin'] = 1

slurm_dict = {}
slurm_dict['email'] = 'eragasa@ufl.edu'
slurm_dict['qos'] = 'phillpot'
slurm_dict['ntasks'] = 16
slurm_dict['time'] = "1:00:00"

# not sure if this works
#poscar = vasp.Poscar(structure_filename)
poscar = vasp.Poscar()
poscar.read(structure_filename)

# this approach is taken from the VASP methodology of generating automated kpoint meshes
# Ref: https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html

linear_kpoint_density = 3 # kpoints/(units reciprocal lattice)

def get_kpoint_mesh(simulation_cell,linear_kpoint_density):
    """
    Args:
        simulation_cell(pypospack.crystal.SimulationCell): self-explanatory. 
    """
    # assert isinstance(simulation_cell,crystal.SimulationCell)
    b1_length = np.dot(simulation_cell.b1,simulation_cell.b1)**0.5
    b2_length = np.dot(simulation_cell.b2,simulation_cell.b2)**0.5
    b3_length = np.dot(simulation_cell.b3,simulation_cell.b3)**0.5
    
    #linear_kpoint_density = kpoints_1/b1_length -->i
    kpoints_1 = linear_kpoint_density*b1_length
    kpoints_2 = linear_kpoint_density*b2_length
    kpoints_3 = linear_kpoint_density*b3_length
    
    kpoints_1_round = int(kpoints_1+0.5)
    kpoints_2_round = int(kpoints_2+0.5)
    kpoints_3_round = int(kpoints_3+0.5)

    return [kpoints_1_round,kpoints_2_round,kpoints_3_round]

kpoint_meshes = {}
def kpoint_mesh_to_str(kpoint_mesh):
    return "kp_{}_{}_{}".format(
        kpoint_mesh[0],
        kpoint_mesh[1],
        kpoint_mesh[2])

old_kpoint_mesh = None
new_kpoint_mesh = None
for linear_kpoint_density in np.arange(1,10,0.1):
    old_kpoint_mesh = copy.copy(new_kpoint_mesh)
    new_kpoint_mesh = get_kpoint_mesh(
        poscar,
        linear_kpoint_density)
    print(linear_kpoint_density,':',
          new_kpoint_mesh)
    if old_kpoint_mesh is not None:
        kpoints_are_different = not all([
            old_kpoint_mesh[0] == new_kpoint_mesh[0],
            old_kpoint_mesh[1] == new_kpoint_mesh[1],
            old_kpoint_mesh[2] == new_kpoint_mesh[2]])
        if kpoints_are_different:
            key = kpoint_mesh_to_str(new_kpoint_mesh)
            kpoint_meshes[key] = list(new_kpoint_mesh)
            print('\tadd')       

