import pytest
import shutil
import pypospack.task.vasp as tsk_vasp

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

cleanup(directory)
converg_encut = tsk_vasp.VaspEncutConvergence(
	directory=directory,
	structure=structure_filename,
	xc=xc,
	incar_dict=incar_dict,
        slurm_dict=slurm_dict,
        full_auto=full_auto)

print("VaspEncutConvergence.directory:\n\t{}".format(
        converg_encut.directory))
print("VaspEncutConvergence.orig_directory:\n\t{}".format(
        converg_encut.orig_directory))
print("VaspEncutConvergence.structure_file:\n\t{}".format(
        converg_encut.structure_filename))
print("VaspEncutConvergence.filename_results:\n\t{}".format(
        converg_encut.filename_results))

converg_encut.create_simulations()

#--------
# full-auto test
cleanup(directory)
full_auto=True
converg_encut = tsk_vasp.VaspEncutConvergence(
	directory=directory,
	structure=structure_filename,
	xc=xc,
	incar_dict=incar_dict,
        slurm_dict=slurm_dict,
        full_auto=full_auto)

