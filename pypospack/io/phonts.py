# -*- coding: utf-8 -*-
"""pypospack.io.phonts

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

Job Schedulers
==============

Due to the wide variety and implementations of the cluster job schedulers,
it is not possible to write all the possibilities of the job scheduler. The
general strategy for dealing with job submissions is to create a 

pypospack.io.slurm
pypospack.io.rocks
pypospack.io.pbs

"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, shutil, pathlib, copy, shelex, subprocess,sys
import numpy as np
import pypospack.io.vasp as vasp
import pypospack.io.slurm as slurm
import pypospack.crystal as crystal
import pypospack.potential as potential
import pypospack.crystal as crystal

# *****************************************************************************
# ****    ERROR EXCEPTION HANDLING CLASSES                                 ****
# *****************************************************************************

class PhontsError(Exception):
    def __init__(self,*args,**kwargs):
        """Error class for reading/writing VASP INCAR IO issues """
        Exception.__init__(self,*args,**kwargs)

# *****************************************************************************
# ****     CORE CLASSES                                                    ****
# *****************************************************************************

phonts_fp_interfaces = ['VASP','QE','LAMMPS']
job_schedulers = ['slurm','rocks','pbs']

def simulation_cell_to_gulp_string(sim_cell):
    if not isinstance(sim_cell,crystal.SimulationCell):
        err_msg = 'simulation cell is not an instance of '
        err_msg += 'pypospack.crystal.SimulationCell'
        raise ValueError(err_msg)

    # check if system is cubic
    a0 = sim_cell.a0
    H = sim_cell.H
    cond1 = H[0,0] != 0 and H[0,1] == 0 and H[0,2] == 0
    cond2 = H[1,0] == 0 and H[1,1] != 1 and H[1,2] == 0
    cond3 = H[2,0] == 0 and H[2,1] == 0 and H[2,2] == 0
    if cond1 and cond2 and cond3:
        pass
    else:
        err_msg = 'simulation cell is not orthogonal'
        raise ValueError(err_msg)

def gulp_phonon_section_to_string(self,shrink=[8,8,8],kpoints=[10,10,10]):
    str_out = (\
        "shrink\n"
        "{} {} {}\n"
        "kpoints\n"
        "{} {} {}\n"
        "output freq text freq.gulp 12 \n"
        "output phonon text phonon.gulp\n"
        "output osc test phonon.osc\n").format(\
                shrink[0],shrink[1],shrink[2],
                kpoints[0],kpoints[1],kpoints[2])

    return str_out

def gulp_positions_to_string(self,structure='POSCAR'):
    """ returns the gulp position section to create simulation cell

    Args:
        structure (str): the filename of the POSCAR format file.  Default 
            is 'POSCAR'.  If an object which subclasses the 
            pyposmat.crystal.SimulationCell class is passed, this method
            will use this that class instead.

    Returns
        str:
    """

    sim_cell = None
    if isinstance(structure,crystal.SimulationCell):
        sim_cell = crystal.SimulationCell(poscar)
    else:
        sim_cell = vasp.Poscar()
        try:
            sim_cell.read(structure)
        except FileNotFoundError as e:
            msg = 'PYPOSPACK_ERROR: cannot find the POSCAR file:{}'.format(structure)
            raise

    H = sim_cell.H * sim_cell.a0
    str_out = "vectors\n"
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[0,0],H[0,1],H[0,2])
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[1,0],H[1,1],H[1,2])
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[2,0],H[2,1],H[2,2])
    str_out += "fractional\n"
    for s in sim_cell.symbols:
        for a in sim_cell.atomic_basis:
            if a.symbol == s:
                try:
                    str_out += "{} core {} {} {}\n".format(\
                            s,a.position[0],a.position[1],a.position[2])
                except:
                    print(s)
                    print(a.symbol)
                    print(a.position)
                    raise
    return str_out

def get_charge(symbol):
    return 0

def phonts_kpoints_section(k1=9,k2=9,k3=9,kbte=1):
    """ Define the kpoint mesh for Phonons

    Defines the k-point mesh on which phonon properties are to be calculated.
    kx, ky, and kz determine resolutions in the reciprocal space.  PhonTS
    automatically create an odd-valued grid, even if an even-grid is provided.

    Args:
        kx(int): kpoint mesh resolution in the b1 direction, reciprocal lattice
        ky(int): kpoint mesh resolution in thee b2 direction, reciprocal lattice
        kz(int): kpoint mesh resolution in the b3 direction, reciprocal lattice
        kbte(int): additional resolution for the finer grid in the z direction,
            for more accurate determination of the energy conservation surface
            in thermal conductivity calculations.
    """

    str_out = "kpoints {} {} {} {}".format(k1,k2,k3,kbte)

def phonts_bte_section():
    pass

def phonts_dispersions():
    pass

def phonts_ab_initio_section(is_start=True,
        numerical_2der=True, numerical_3dr=True,
        fp_interface='VASP',asr_3der=False,d3_cutoff=10.0):
    """
    Args:
        numerical_2der(bool):
        numerical_3der(bool):
        d3_cutoff(float):
        fp_interface(str): first-principals interface.  Can be VASP, QE or
            LAMMPS.  Default is set to LAMMPS.
    """

    str_out = ''
    if fp_interface not in phonts_fp_interfaces:
        if fp_interface is not None:
            msg_out = "fp_interface must either be VASP, QE, or LAMMPS or None"
            raise ValueError(msg_out)
    else:
        str_out += 'FP_interface {}'.format(fp_interface)

    if is_start:
        str_out += 'AbInitio T F\n'
    else:
        str_out += 'AbInitio F T\n'



    if numerical_2der is True:

        str_out += 'numerical_2der T\n'
    else:
        str_out += 'numerical_2der F\n'

    if numerical_3der is True:
        str_out += 'numerical_3der T\n'
    else:
        str_out += 'numerical_3der F\n'

    str_out += 'D3_cutoff {}\n'.format(d3_cutoff)

    return str_out

# *****************************************************************************
# ****     CORE CLASSES                                                    ****
# *****************************************************************************
class PhontsSimulation(object):
    """
    Args:
        task_name(str)
        task_dir(str)

    Attributes:
        task_name(str)
        task_dir(str)
        structure_filename(str): Attribute initialized to None.
        structure(pypospack.crystal.SimulationCell): Attribute initialized to 
            None
        simulation_cell(pypospack.io.vasp.Poscar)
        filename(str): name of the phonts input file  Initialized to be
            <task_dir>/phonons_input.dat

        potential(pypospack.potential.Potential)
        fp_interface(str): must be
            VASP, QE, or LAMMPS
    """
    def __init__(self,task_name='phonts',task_dir=None):

        if isinstance(task_name,str):
            self.task_name = task_name
        else:
            err_msg = "task_name must be a string"
            raise ValueError(err_msg)

        if task_dir is None:
            self.task_directory = os.path.join(
                    os.getcwd(),'phonts')
        elif isinstance(task_dir,str):
            self.task_directory = task_dir
        else:
            err_msg = "task_dir must be a string"
            raise ValueError(err_msg)

        self.structure_filename = None
        self.structure = None
        self.simulation_cell = vasp.Poscar()
        self.filename = os.path.join(
                self.task_directory,
                'phonons_input.dat')

        # internal interatomic potential parameters
        self.potential = None
        self.phonts_potential_type = None
        self.phonts_potential_params = None

        # external force evaluation
        self.fp_interface = None
        self.phonts_bin = os.environ['PHONTS_BIN']

        # job manager attributes
        self.job_scheduler = None
        self.slurm_phonts_dict = None
        self.slurm_vasp_dict = None
        self._create_task_directory()

    def write_input_file(self):
        str_out = "# PhonTS input file created pypospack.\n"
        str_out += "# {:*^78}\n".format('force evaluation')
        str_out += self.force_evaluation_to_string()
        str_out += "# {:*^78}\n".format('simulation cell')
        str_out += self.simulation_cell_to_string()
        str_out += "\nend\n"
           
        with open(self.filename,'w') as f:
            f.write(str_out)
        return str_out
    
    def write_submission_script(self):

        if self.job_scheduler == 'slurm':
            self.slurm_phonts_dict['filename'] = os.path.join(
                    self.task_directory,
                    self.slurm_phonts_dict['filename'])
            slurm.write_phonts_batch_script(\
                    filename=self.slurm_phonts_dict['filename'],
                    job_name=self.slurm_phonts_dict['job_name'],
                    email=self.slurm_phonts_dict['email'],
                    qos=self.slurm_phonts_dict['qos'],
                    ntasks=self.slurm_phonts_dict['ntasks'],
                    time=self.slurm_phonts_dict['time'])
        elif self.job_scheduler == 'sge':
            raise NotImplementedError()
        elif self.job_scheduler == 'pbs':
            raise NotImplementedError()

    def run(self):
        self.orig_dir = os.getcwd()
        os.chdir(self.task_directory)
        if self.job_scheduler == 'slurm':
            cmd_s = 'sbatch {}'.format(self.slurm_phonts_dict['filename'])
            args = shlex.split(cmd_s)
            p = subprocess.Popen(args)
    def force_evaluation_to_string(self):
        str_out = ''
        if  self.fp_interface in phonts_fp_interfaces:
            str_out += self._get_external_force_evaluation_to_string() 
        elif self.potential is not None:
            str_out += self._get_internal_force_evaluation_to_string()
        else:
            err_msg = "either the potential must be specified or an "
            err_msg += "external method must be provided"
            raise ValueError(err_msg)
        return str_out
    
    def _get_external_force_evaluation_to_string(self):
        str_out = ''
        if self.fp_interface == 'VASP':
            pass
        elif self.fp_interface == 'QE':
            raise NotImplementedError('fp_interface == QE not supported')
        elif self.fp_interface == 'LAMMPS':
            raise NotImplementedError('fp_interface == LAMMPS not supported')
        return str_out

    def _get_internal_force_evaluation_to_string(self):
        str_out = ''
        str_out += ' '.join([str(v) for v in self.phonts_potential_type]) + '\n'
        str_out += ' '.join([str(v) for v in self.phonts_potential_params]) + '\n'
        return str_out 
    
    def simulation_cell_to_string(self):
        if isinstance(self.structure_filename,str):
            self.structure = vasp.Poscar()
            self.structure.read(self.structure_filename)

        str_out = ""
        if isinstance(self.structure,crystal.SimulationCell):
            # get relevant information
            n_species = len(self.structure.symbols)
            a0 = self.structure.a0
            a1 = self.structure.a1/a0
            a2 = self.structure.a2/a0
            a3 = self.structure.a3/a0
            n_atoms = self.structure.n_atoms

            # build string
            str_out += "species {}\n".format(n_species)
            for s in self.structure.symbols:
                s_charge = self._get_charge(s)
                s_amu = crystal.get_amu(s)
                str_out += "{} {} {}\n".format(s,s_charge,s_amu)
            str_out += "Lattice {:.10f}\n".format(a0)
            str_out += "cell {:.10f} {:.10f} {:.10f}\n".format(a1,a2,a3)
            str_out += "natoms {}\n".format(n_atoms)
            for a in self.structure.atomic_basis:
                s = a.symbol
                x = a.position[0]
                y = a.position[1]
                z = a.position[2]
                str_out += "{} {:10f} {:10f} {:10f}\n".format(s,x,y,z)
        else:
            msg_out = "structure must be an instance of "
            msg_out += "pypospack.crystal.SimulationCell"
            raise ValueError(msg_out)

        return str_out
    def _get_charge(self,s):
        return 0
    def _create_task_directory(self):
        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)
        os.mkdir(self.task_directory)


