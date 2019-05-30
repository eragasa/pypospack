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
            
if __name__ == '__main__':
    
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
