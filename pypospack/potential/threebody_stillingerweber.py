# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from collections import OrderedDict
from pypospack.potential import Potential

class StillingerWeberPotential(Potential):
    def __init__(self,symbols):
        """
        Args:
            symbols: list of string
        Attributes:
            symbols
            potential_type
            is_charge
        References:
            http://lammps.sandia.gov/doc/pair_sw.html
        """
        _potential_type = 'stillingerweber'
        _is_charge = False

        Potential.__init__(self,
                symbols=symbols,
                potential_type=_potential_type,
                is_charge=_is_charge)

        self.lmps_parameter_filename = "lmps_parameter_filename"

    def _init_parameter_names(self):
        # TODO: This is only written for a single element potential
        _symbols = self.symbols
        _n_symbols = len(_symbols)
        for i in range(_n_symbols):
            for j in range(_n_symbols):
                for k in range(_n_symbols):
                    el1 = _symbols[i]
                    el2 = _symbols[j]
                    el3 = _symbols[k]
                    self._add_parameter_names(el1,el2,el3)

    def _add_parameter_names(self,el1,el2,el3):
        s = "{}{}{}".format(el1,el2,el3)
        
        sw_param_names = [
                'epsilon',
                'sigma',
                'a',
                'lambda',
                'gamma',
                'costheta0',
                'A',
                'B',
                'p',
                'q',
                'tol'
            ]
        
        self.parameter_names = []
        for p in sw_param_names:
            self.parameter_names.append("{}_{}".format(s,p))

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def lammps_potential_section_to_string(self,parameters=None):

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        fname_params= self.lmps_parameter_filename

        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        str_out += "pair_style sw\n"
        for s in self.symbols:
            str_out += "pair_coeff * * {} {}\n".format(fname_params,s)
        str_out += "\n"

        #<--------- neighbor lists moved to pypospack.task.lammps.LammpsTask
        #str_out += "neighbor 1.0 bin\n"
        #str_out += "neigh_modify every 1 delay 0 check yes\n"

        return str_out
   
    def write_lammps_parameter_file(self,dst_dir,dst_filename):
        assert isinstance(dst_dir,str)
        assert isinstance(dst_filename,str)

        _strout = self.lammps_paramfile_file_to_string()
        with open(os.path.join(dst_dir,dst_filename)) as f:
            f.write(_strout)

    def lammps_parameter_file_to_string(self,parameters=None):
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]
        
        str_out = ''
        for i, s1 in enumerate(self.symbols):
            for j, s2 in enumerate(self.symbols):
                for k, s3 in enumerate(self.symbols):
                    s = '{}{}{}'.format(s1,s2,s3)
                    str_out += '{} {} {}'.format(s1,s2,s3)
                    str_out += ' ' + str(self.parameters['{}_epsilon'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_sigma'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_a'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_lambda'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_gamma'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_costheta0'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_A'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_B'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_p'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_q'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_tol'.format(s)])
                    str_out += '\n'

        return str_out
       

