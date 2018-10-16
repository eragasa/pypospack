import os
from collections import OrderedDict
import pypospack.potential as potential

from pypospack.task.lammps import LammpsSimulationError, LammpsSimulation

class LammpsElasticCalculation(LammpsSimulation):
    """ Class for LAMMPS elastic calculation

    This class sets up, runs and processes the relaxation of atomic poistions
    to find the lowest energy struction in the local basin of attraction. The
    simulation cell is then deformed and atomic positions relaxed to calculate
    the energy differences to calculate elements of the elastic tensor

    This class is based on the LAMMPS script written by Aiden Thompson

    Args:
        task_name (str): the name of the task
        task_directory (str): the directory of the task
    """
    def __init__(self,
            task_name,
            task_directory,
            structure_filename,
            restart=False,
            fullauto=False):

        _task_type = 'lmps_elastic'
        LammpsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                task_type=_task_type,
                structure_filename=structure_filename,
                restart=restart,
                fullauto=fullauto)
        
        self.lammps_results_names = [
                'c11','c12','c13',
                'c22','c23','c24','c25','c26',
                'c33','c34','c35','c36',
                'c44','c45','c46',
                'c55','c56',
                'c66']

    def postprocess(self):
        LammpsSimulation.postprocess(self)

    def on_post(self,configuration=None):
        #self.__get_results_from_lammps_outputfile()
        self.get_elastic_tensor()
        LammpsSimulation.on_post(self,configuration=configuration)

    def on_finished(self,configuration=None): 
        pass

    def write_potential_file(self):
        if self.potential is None:
            return
        
        _setfl_dst_filename = None
        if isinstance(self.potential,potential.EamPotential):
            _setfl_dst_filename = os.path.join(
                    self.task_directory,
                    "{}.eam.alloy".format("".join(self.potential.symbols)))
            _str_out = self.potential.lammps_potential_section_to_string(
                setfl_dst_filename=_setfl_dst_filename)
        
        # <-------- FOR STILLINGER WEBER POTENTIALS
        elif isinstance(self.potential,potential.StillingerWeberPotential):
            # define the filename --- SiO.parameters, Si.parameters
            _symbols_str = "".join(self.potential.symbols)
            _p_fname = "{}.parameters".format(_symbols_str)

            # set the name of the output file
            self.potential.lmps_parameter_filename = _p_fname
            
            # get the string of potential.mod
            _str_out = self.potential.lammps_potential_section_to_string()

            # write out the potential parameter file
            _str_lmps_params = self.potential.lammps_parameter_file_to_string()
            
            _p_fname_dst = os.path.join(self.task_directory,_p_fname)
            with open(_p_fname_dst,'w') as f:
                f.write(_str_lmps_params)

        #default behavior
        else:
            _str_out = self.potential.lammps_potential_section_to_string()

        _str_out += "\n"
        
        # coulumbic charge summation         
        if self.potential.is_charge:
            _str_out += "kspace_style pppm 1.0e-5\n"
            _str_out += "\n"
        
        # neighborlists
        _str_out += "neighbor 1.0 bin\n"
        _str_out += "neigh_modify every 1 delay 0 check yes\n"

        # Setup output
        _str_out += "\n".join([
                "thermo	1",
                "thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol",
                "thermo_modify norm no\n"])
        
        _lammps_potentialmod_filename = os.path.join(
                self.task_directory,
                self.lammps_potentialmod_filename)
        
        with open(_lammps_potentialmod_filename,'w') as f:
            f.write(_str_out)

    def write_lammps_input_file(self):
        """ writes LAMMPS input file

        This method is modified from the LammpsSimulation template due to 
        the need for the multiple files.
        
        Args:
            filename (str): name of the input file for LAMMPS. Default is
                'lammps.in'.
        Attributes:
        
        """

        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.task_directory,
                self.lammps_input_filename)
        with open(filename,'w') as f:
            f.write(str_out)

        str_out = self.lammps_init_mod_to_string()
        filename = os.path.join(self.task_directory,'init.mod')
        with open(filename,'w') as f:
            f.write(str_out)

        str_out = self.lammps_displace_mod_to_string()
        filename = os.path.join(self.task_directory,'displace.mod')
        with open(filename,'w') as f:
            f.write(str_out)


    def lammps_input_file_to_string(self):
        str_out = (
            "include init.mod\n"
            "include potential.mod\n"
            "# ---- Compute initial state\n"
            "fix 3 all box/relax aniso 0.0\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "variable tmp equal pxx\n"
            "variable pxx0 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy0 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz0 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz0 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz0 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy0 equal ${tmp}\n"
            "\n"
            "variable tmp equal lx\n"
            "variable lx0 equal ${tmp}\n"
            "variable tmp equal ly\n"
            "variable ly0 equal ${tmp}\n"
            "variable tmp equal lz\n"
            "variable lz0 equal ${tmp}\n"
            "\n"
            " # ---- define the derivatives w.r.t. strain components\n"
            "variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}\n"
            "variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}\n"
            "variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}\n"
            "\n"
            "# ---- write restart files\n"
            "unfix 3\n"
            "write_restart restart.equil\n"
            "# ---- uxx Perturbation\n"
            "variable dir equal 1\n"
            "include displace.mod\n"
            "# ---- uyy Perturbation\n"
            "variable dir equal 2\n"
            "include displace.mod\n"
            "# ---- uzz Perturbation\n"
            "variable dir equal 3\n"
            "include displace.mod\n"
            "# ---- uyz Perturbation\n"
            "variable dir equal 4\n"
            "include displace.mod\n"
            "# ---- uxz Perturbation\n"
            "variable dir equal 5\n"
            "include displace.mod\n"
            "# ---- uxy Perturbation\n"
            "variable dir equal 6\n"
            "include displace.mod\n"
            "\n"
            "# ---- Output final values\n"
            "variable C11all equal ${C11}\n"
            "variable C22all equal ${C22}\n"
            "variable C33all equal ${C33}\n"
            "variable C12all equal 0.5*(${C12}+${C21})\n"
            "variable C13all equal 0.5*(${C13}+${C31})\n"
            "variable C23all equal 0.5*(${C23}+${C32})\n"
            "variable C44all equal ${C44}\n"
            "variable C55all equal ${C55}\n"
            "variable C66all equal ${C66}\n"
            "variable C14all equal 0.5*(${C14}+${C41})\n"
            "variable C15all equal 0.5*(${C15}+${C51})\n"
            "variable C16all equal 0.5*(${C16}+${C61})\n"
            "variable C24all equal 0.5*(${C24}+${C42})\n"
            "variable C25all equal 0.5*(${C25}+${C52})\n"
            "variable C26all equal 0.5*(${C26}+${C62})\n"
            "variable C34all equal 0.5*(${C34}+${C43})\n"
            "variable C35all equal 0.5*(${C35}+${C53})\n"
            "variable C36all equal 0.5*(${C36}+${C63})\n"
            "variable C45all equal 0.5*(${C45}+${C54})\n"
            "variable C46all equal 0.5*(${C46}+${C64})\n"
            "variable C56all equal 0.5*(${C56}+${C65})\n"
            "\n"
            "print \"c11 = ${C11all} ${cunits}\"\n"
            "print \"c22 = ${C22all} ${cunits}\"\n"
            "print \"c33 = ${C33all} ${cunits}\"\n"
            "print \"c12 = ${C12all} ${cunits}\"\n"
            "print \"c13 = ${C13all} ${cunits}\"\n"
            "print \"c23 = ${C23all} ${cunits}\"\n"
            "print \"c44 = ${C44all} ${cunits}\"\n"
            "print \"c55 = ${C55all} ${cunits}\"\n"
            "print \"c66 = ${C66all} ${cunits}\"\n"
            "print \"c14 = ${C14all} ${cunits}\"\n"
            "print \"c15 = ${C15all} ${cunits}\"\n"
            "print \"c16 = ${C16all} ${cunits}\"\n"
            "print \"c24 = ${C24all} ${cunits}\"\n"
            "print \"c25 = ${C25all} ${cunits}\"\n"
            "print \"c26 = ${C26all} ${cunits}\"\n"
            "print \"c34 = ${C34all} ${cunits}\"\n"
            "print \"c35 = ${C35all} ${cunits}\"\n"
            "print \"c36 = ${C36all} ${cunits}\"\n"
            "print \"c45 = ${C45all} ${cunits}\"\n"
            "print \"c46 = ${C46all} ${cunits}\"\n"
            "print \"c56 = ${C56all} ${cunits}\"\n"
            "print \"lammps_sim_done\"\n")
        return str_out

    def lammps_init_mod_to_string(self):
        if self.potential.is_charge:
            _atom_style = 'charge'
        else:
            _atom_style = 'atomic'
        _structure_file = self.lammps_structure_filename
        
        str_out = (\
                "# ---- init.mod file\n"
                "variable up equal 1.0e-6\n"
                "units metal\n"
                "dimension 3\n"
                "boundary p p p\n"
                "atom_style {atom_style}\n"
                "atom_modify map array\n"
                "variable cfac equal 1.0e-4\n"
                "variable cunits string GPa\n"
                "# ---- define minimization parameters\n"
                "variable etol equal 0.0\n"
                "variable ftol equal 1.0e-10\n"
                "variable maxiter equal 100\n"
                "variable maxeval equal 1000\n"
                "variable dmax equal 1.0e-2\n"
                "# --- read data structure\n"
                "read_data {structure_file}\n"
            ).format(\
                atom_style=_atom_style,
                structure_file=_structure_file
            )
        return str_out

    def lammps_displace_mod_to_string(self):
        str_out = (\
            "# NOTE: This script should not need to be\n"
            "# modified. See in.elastic for more info.\n"
            "# Find which reference length to use\n"
            "\n"
            "if \"${dir} == 1\" then &\n"
            "   \"variable len0 equal ${lx0}\"\n" 
            "if \"${dir} == 2\" then &\n"
            "   \"variable len0 equal ${ly0}\"\n" 
            "if \"${dir} == 3\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 4\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 5\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 6\" then &\n"
            "   \"variable len0 equal ${ly0}\"\n" 
            "\n"
            "# Reset box and simulation parameters\n"
            "\n"
            "clear\n"
            "read_restart restart.equil remap\n"
            "include potential.mod\n"
            "\n"
            "# Negative deformation\n"
            "\n"
            "variable delta equal -${up}*${len0}\n"
            "if \"${dir} == 1\" then &\n"
            "   \"change_box all x delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 2\" then &\n"
            "   \"change_box all y delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 3\" then &\n"
            "   \"change_box all z delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 4\" then &\n"
            "   \"change_box all yz delta ${delta} remap units box\"\n"
            "if \"${dir} == 5\" then &\n"
            "   \"change_box all xz delta ${delta} remap units box\"\n"
            "if \"${dir} == 6\" then &\n"
            "   \"change_box all xy delta ${delta} remap units box\"\n"
            "\n"
            "# Relax atoms positions\n"
            "\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "# Obtain new stress tensor\n"
            "\n"
            "variable tmp equal pxx\n"
            "variable pxx1 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy1 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz1 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy1 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz1 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz1 equal ${tmp}\n"
            "\n"
            "# Compute elastic constant from pressure tensor\n"
            "\n"
            "variable C1neg equal ${d1}\n"
            "variable C2neg equal ${d2}\n"
            "variable C3neg equal ${d3}\n"
            "variable C4neg equal ${d4}\n"
            "variable C5neg equal ${d5}\n"
            "variable C6neg equal ${d6}\n"
            "\n"
            "# Reset box and simulation parameters\n"
            "\n"
            "clear\n"
            "read_restart restart.equil remap\n"
            "include potential.mod\n"
            "\n"
            "# Positive deformation\n"
            "\n"
            "variable delta equal ${up}*${len0}\n"
            "if \"${dir} == 1\" then &\n"
            "   \"change_box all x delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 2\" then &\n"
            "   \"change_box all y delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 3\" then &\n"
            "   \"change_box all z delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 4\" then &\n"
            "   \"change_box all yz delta ${delta} remap units box\"\n"
            "if \"${dir} == 5\" then &\n"
            "   \"change_box all xz delta ${delta} remap units box\"\n"
            "if \"${dir} == 6\" then &\n"
            "   \"change_box all xy delta ${delta} remap units box\"\n"
            "\n"
            "# Relax atoms positions\n"
            "\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "# Obtain new stress tensor\n"
            "\n"
            "variable tmp equal pe\n"
            "variable e1 equal ${tmp}\n"
            "variable tmp equal press\n"
            "variable p1 equal ${tmp}\n"
            "variable tmp equal pxx\n"
            "variable pxx1 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy1 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz1 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy1 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz1 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz1 equal ${tmp}\n"
            "\n"
            "# Compute elastic constant from pressure tensor\n"
            "\n"
            "variable C1pos equal ${d1}\n"
            "variable C2pos equal ${d2}\n"
            "variable C3pos equal ${d3}\n"
            "variable C4pos equal ${d4}\n"
            "variable C5pos equal ${d5}\n"
            "variable C6pos equal ${d6}\n"
            "\n"
            "# Combine positive and negative\n"
            "\n"
            "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})\n"
            "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})\n"
            "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})\n"
            "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})\n"
            "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})\n"
            "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})\n"
            "\n"
            "# Delete dir to make sure it is not reused\n"
            "\n"
            "variable dir delete\n")
        return str_out

    def get_elastic_tensor(self):
        filename = os.path.join(self.task_directory,'lammps.out')
        with open(filename,'r') as f:
            lines = f.readlines()

        _lammps_results_names = self.lammps_results_names

        self.results = OrderedDict()


        _task_name = self.task_name
        for i,line in enumerate(lines):
            for name in _lammps_results_names:
                if line.startswith('{} = '.format(name)):
                    # get the result name and result value from the file
                    # _rk = key to store result
                    # _rv = value to store
                    _rk = '{}.{}'.format(_task_name,name)
                    _rv = float(line.split('=')[1].split()[0].strip())

                    self.results[_rk] = _rv
                elif line.startswith('ERROR'):
                    print('name:{}'.format(name))
                    print('line:{}'.format(line.strip))
                    raise NotImplementedError
        
        #_all_components_exist = all([n in self.results for n in _lammps_results_names])
        #if not _all_components_exist:
        #    for n in _lammps_results_names:
        #        print(n,n in self.results.keys())


