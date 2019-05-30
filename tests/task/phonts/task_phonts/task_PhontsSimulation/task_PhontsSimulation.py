"""
pypospack.io.phonts

This module provides input and output routines and classes for managing
simulations for Phonon Transport Simulation (PhonTS), written by 
A. Chernatiynskiy and S.R. Phillpot.  PhonTS is a lattice dynamics code 
that calculates thermal conductivity via the solution of the Boltzmann 
Transport Equation (BTE) for phonons.

This module's purpose is to automate the work required to interface PhonTS to
other codes such as LAMMPS for additional potentials, and first principals
electronic structure calculations (VASP).  While an interface is defined
for Quantum Expresso QE) for PhonTS, this interface is not implemented with
pypospack.io.phonts.

For a more a more in-depth understand of the PhonTS refer to the foundation
papers:

A. Chernatiynskiy and S.R. Phillpot, PRB 82, 134301 (2010)
A. Chernatiynskiy and S.R. Phillpot, Comp. Phys. Comm. (2015)

This module does not implement all the functionality which is available in
PhonTS.

Requirements
============

The use of this module requires setting the enviornment variable the 
.bash_profile or similar script on your system

export PHONTS_BIN=~/bin/PhonTS

"""
import sys
import os, copy, shutil, subprocess
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

import pypospack.task.Task

#<<<<<< These classes will end up in pypospack.io.phonts
class PhontsInputFile(object):
    pass

class PhontsOutputFile(object):
    pass



#<<<<<< These classes will end up in pypospack.task.phonts
class PhontsSimulation(pypospack.task.Task):
    def __init__(task_name,task_directory):
        pypospack.task.Task.__init__(task_name,task_directory)

    def set_structure(self,structure):
        # does task.Task do this already for me?
class PhontsDispersion(PhontsSimulation):
    def __init__(task_name,task_directory):
        PhontSimulation.task.Task.__init__(task_name,task_directory)

    def set_structure(self,structure):
        PhontsSimulation.set_structure(structure)
    
    def set_dispersion_path(self,BZ_points,BZ_paths,N_BZ_ppp):
        """
        Args:
            BZ_points(diction of str,float): list of unique high-symmetry point in hte BZ zone
            BZ_paths(list of list of str): list of paths
            N_BZ_ppp(int): numbe of points per path.
        """
        phonts_dispersion_path_cmd = "do_dispersion {n_paths} {n_pts_per_path}"

        all_paths = []
        for bz_path in BZ_paths:
            for i in range(1,len(bz_path)):
                bzp1 = bz_path[i-1]
                bzp2 = bz_path[i]
                if all([ [bzp1,bzp2] not in all_paths,
                         [bzp2,bzp1] not in all_paths]):
                    all_paths.append([bzp1,bzp2])
        
        n_paths = len(all_paths)

        str_out = "do_dispersion {} {}".format(
                n_paths,
                N_BZ_points_per_path)
        for path in in_paths:
            "{
if __name__ == '__main__':
    
    # this was generated from www.aflowlib.org/aflow_online.html
    BZ_points = OrderedDict()
    BZ_points['G']=[0.000,0.000,0.000]
    BZ_points['X']=[0.500,0.000,0.500]
    BZ_points['W']=[0.500,0.250,0.750]
    BZ_points['K']=[0.375,0.375,0.750]
    BZ_points['L']=[0.500,0.500,0.500]
    BZ_points['U']=[0.625,0.250,0.625]
    BZ_paths = []
    BZ_paths.append(['G','X','W','K','G'])
    BZ_paths.append(['G','L','U','W','L','K'])
    BZ_paths.append(['U','X'])
    BZ_points_per_path
    vasp_filename = 'Si_dia.gga.relax.vasp'
    vasp_input_filename = os.path.join(os.getcwd(),'rsrc',vasp_filename)

    #### TEST INITIALIZE ####
    task_name = 'phonts'
    task_directory = os.path.join(\
            os.getcwd(),
            task_name)
    task = PhontsPhononCalculation(task_name,task_directory)
    assert task.task_name == task_name
    assert task.task_directory == task_directory

    #### get_atomic_structure_from_filename ####
    sim_cell = vasp.Poscar()
    sim_cell.read(vasp_input_filename)
    str_out = "# --- CRYSTAL STRUCTURE ---\n"
    str_out += "species {}".format(len(sim_cell.symbols))
    for s in sim_cell.symbols:
        s_charge = get_charge(s)
        s_amu = crystal.get_amu(s)
        str_out += "{} {} {}\n"(s, s_charge, s_amu)
    str_out += "Lattice {:.10f}\n".format(sim_cell.a0)
    str_out += "cell {:10f} {:10f} {:10f}\n".format(sim_cell.a1,sim_cell.a2,sim_cell.a3)
    str_out += "natoms {}\n".format(sim_cell.n_atoms)
    for a in sim_cell.atomic_basis:
        str_out += "{} {:10f} {:10f} {:10f}\n".format(
                a.symbol, a.position[0], a.position[1], a.position[2])
    print(str_out)
    sys.exit()

    #### TEST POTENTIAL POTENTIAL ####
    print('----- test that buckingham potential provides the right format -----')
    task.potential = potential.Buckingham(['Mg','O'])
    task.param_dict = copy.deepcopy(param_dict)
    print(task.potential.gulp_potential_section_to_string(param_dict))

    #### TEST IF WE CAN CAN WRITE THE INPUT FILE ####
    gulp_input_filename = os.path.join(task.task_directory,'gulp.in')
    poscar_input_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    task.write_gulp_input_file(\
            filename=gulp_input_filename,
            poscar=vasp_input_filename)

    #### TEST IF WE CAN RUN THE BUCKINGHAM POTENTIAL ####
    task.structure_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    task.run()


    task = GulpPhononCalculation(task_name,task_directory)
    task.potential = potential.Buckingham(['Mg','O'])
    task.param_dict = copy.deepcopy(param_dict)
    task.write_gulp_input_file(\
            filename=gulp_input_filename,
            poscar=vasp_input_filename)
    task.run()
