
"""
This script tests the functionality of the pypospack.io.phonts module by 
recreating the example in the PhonTS example/argon

"""
import copy
import numpy as np
import pypospack.io.phonts as phonts
import pypospack.crystal as crystal

# define simulation cell directly
ar_fcc = crystal.SimulationCell()
ar_fcc.a0 = 3.9620
ar_fcc.add_atom('Ar',[0.0,0.0,0.0])
ar_fcc.add_atom('Ar',[0.5,0.5,0.0])
ar_fcc.add_atom('Ar',[0.5,0.0,0.5])
ar_fcc.add_atom('Ar',[0.0,0.5,0.5])

# define simulation cell through poscar file
# poscar_filename = 'POSCAR'
# acc_fcc = vasp.Poscar()
# acc_fcc.read('POSCAR')

# this is a leonnard jones potential
# TODO: LJ is not currently implemented in pypospack.potentials
# TODO: pypospack.potentials.LeonnardJones.to_phonts_string() not implemented
phonts_potential_type = ['exp-6',1]
phonts_potential_params = ['Ar','Ar',0.010531316,13.,3.85,5.625,9.625]

slurm_phonts_dict = {}
slurm_phonts_dict['filename'] = 'slurm_submit.sh'
slurm_phonts_dict['job_name'] = 'argon_test'
slurm_phonts_dict['email'] = 'eragasa@ufl.edu'
slurm_phonts_dict['qos'] = 'phillpot'
slurm_phonts_dict['ntasks'] = 16
slurm_phonts_dict['time'] = '1:00:00'

if __name__ == '__main__':
    phonts_ar_fcc = phonts.PhontsSimulation()
    phonts_ar_fcc.structure = ar_fcc

    # set an internal potential
    phonts_ar_fcc.potential = "LJ"
    phonts_ar_fcc.phonts_potential_type = phonts_potential_type
    phonts_ar_fcc.phonts_potential_params = phonts_potential_params

    # test creating a simulation cell to string()
    # print('# {:*^78}'.format('SIMULATION_CELL'))
    # structure_str = phonts_ar_fcc.simulation_cell_to_string()
    # print(structure_str)

    # create input file()
    phonts_ar_fcc.write_input_file()
    # testing code below
    # print('# {:*^78}'.format('SIMULATION_CELL'))
    # input_file_str = phonts_ar_fcc.write_input_file()
    # print(input_file_str)

    # create submission script()
    phonts_ar_fcc.job_scheduler = 'slurm'
    phonts_ar_fcc.slurm_phonts_dict = copy.deepcopy(slurm_phonts_dict)
    phonts_ar_fcc.write_submission_script()

