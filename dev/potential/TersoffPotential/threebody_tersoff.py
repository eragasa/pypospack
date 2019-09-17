# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os
from collections import OrderedDict
from pypospack.potential import ThreeBodyPotential


def get_3body_parameter_names(symbols, names):

    parameter_names = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            if i1 <= i2:
                for n in symbols:
                    parameter_names.append('{}{}.{}'.format(s1,s2,n))
            for s3 in symbols:
                for n in names:
                    parameter_names.append('{}{}{}.{}'.format(s1,s2,s3,n))
    return parameter_names

class TersoffPotential(ThreeBodyPotential):

    potential_type = 'tersoff'
    

    twobody_parameter_names = ['n', 'beta', 'lambda2', 'B', 'lambda1', 'A', 'R', 'D']
    threebody_parameter_names= ['m', 'gamma', 'lambda3', 'c', 'd', 'costheta0']

    def __init__(self, symbols):
        """
        Args:
            symbols: list of str
        Attributes:
            symbols
            potential_type
            is_charge
        References:
        """

        potential_type = TersoffPotential.potential_type
        is_charge = False

        ThreeBodyPotential.__init__(self,
                                    symbols=symbols,
                                    potential_type=potential_type,
                                    is_charge=is_charge)
        self.lmps_parameter_filename = 'lmps_parameter_filename'

    def _init_parameter_names(self):

        symbols = self.symbols
        self.parameter_names = []
        for i1,s1 in enumerate(symbols):
            for i2,s2 in enumerate(symbols):
                if i1 <= i2:
                    for n in TersoffPotential.twobody_parameter_names:
                        self.parameter_names.append('{}{}.{}'.format(s1,s2,n))
                for s3 in symbols:
                    for n in TersoffPotential.threebody_parameter_names:
                        self.parameter_names.append('{}{}{}.{}'.format(s1,s2,s3,n))
        return self.parameter_names

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def get_parameter_names(self, symbols=None):
        if symbols is not None:
            self.symbols = symbols
        symbols_ = self.symbols

        parameter_names = []
        for i1,s1 in enumerate(symbols):
            for i2,s2 in enumerate(symbols):
                if i1 <= i2:
                    for n in TersoffPotential.twobody_parameter_names:
                        parameter_names.append('{}{}.{}'.format(s1,s2,n))
                for s3 in symbols:
                    for n in TersoffPotential.threebody:
                        parameter_names.append('{}{}{}.{}'.format(s1,s2,s3,n))

        self.parameter_names = parameter_names
        return parameter_names

    def lammps_potential_section_to_string(self, parameters=None):

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]


        str_out = ''
        for i, s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1, self._get_mass(s))
        str_out += "\n"

        for i, s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s, i+1)
        str_out += "\n"

        parameter_filename_ = self.lmps_parameter_filename
        str_symbols_ = " ".join(self.symbols)
        str_out += "pair_style tersoff\n"
        str_out += "pair_coeff * * {} {}\n".format(parameter_filename_, str_symbols_)
        return str_out

    def write_lammps_parameter_file(self,dst_dir,dst_filename):
        assert type(dst_dir) == str
        assert type(dst_filename) == str
        _strout = self.lammps_parameter_file_to_str()
        with open(os.path.join(dst_dir,dst_filename)) as f:
            f.write(_strout)

    def lammps_parameter_file_to_str(self):
        str_out = ''
        for i, el1 in enumerate(self.symbols):
            for j, el2 in enumerate(self.symbols):
                for k, el3 in enumerate(self.symbols):
                    line_args_names = [
                            'element1', 'element2', 'element3', 
                            'm', 'gamma', 'lambda3', 'c', 'd', 'costheta0', 
                            'n', 'beta', 'lambda2', 'B', 'R', 'D', 'lambda1', 'A']

                    s = '{}{}{}'.format(el1, el2, el3)
                    for n in line_arg_names:
                        if n in TersoffPotential.twobody_parameter_names:
                            pn1 = '{}{}_{}'.format(el1, el2, n)
                            pn2 = '{}{}_{}'.format(el2, el1, n)
                            if pn1 in self.parameters:
                                str_out += ' ' + self.parameters[pn1]
                            elif pn2 in self.parameters:
                                str_out += ' ' + self.parameters[pn2]
                            else:
                                raise NameError()
                        elif n in TersoffPotential.threebody_parameter_names:
                            pn = '{}{}{}_{}'.format(el1, el2, el3, n)
                            str_out += ' ' + self.parameters[pn]
                        else:
                            msg = "cannot find a parameter value for {}{}{}_{}".format(
                                    el1, el2, el3, n
                            )
                            raise ValueError(msg)
                    str_out += '\n'
        return str_out

    def write_lammps_parameter_file(self,dst_dir,dst_filename):

        str_out = self.lammps_parameter_file_to_str()
        filename_ = os.path.join(dst_dir, dst_filename)

        with open(filename_, 'w') as f:
            f.write(str_out)

    def read_lammps_potential_file(self, filename):
        if self.symbols is None:
            self.symbols = []

        parameters = OrderedDict()
        
        with open(filename, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith('#') or line == '\n': 
                pass
            else:
                line_args_names = [
                        'element1', 'element2', 'element3', 
                        'm', 'gamma', 'lambda3', 'c', 'd', 'costheta0', 
                        'n', 'beta', 'lambda2', 'B', 'R', 'D', 'lambda1', 'A']
                line_args = [k.strip() for k in line.strip().split()]
                symbol1 = line_args[0]
                symbol2 = line_args[1]
                symbol3 = line_args[2]
                if symbol1 not in self.symbols: self.symbols.append(symbol1)
                if symbol2 not in self.symbols: self.symbols.append(symbol2)
                if symbol3 not in self.symbols: self.symbols.append(symbol3)
     
                for i,v in enumerate(zip(line_args_names, line_args)):
                    if i > 2:
                        parameter_name = '{}{}{}_{}'.format(symbol1,symbol2,symbol3,v[0])
                        try:
                            parameter_value = int(v[1])
                        except ValueError as e:
                            parameter_value = float(v[1])
                        parameters[parameter_name] = parameter_value

        self.parameters = OrderedDict()
        for i1, s1 in enumerate(self.symbols):
            for i2, s2 in enumerate(self.symbols):
                for i3, s3 in enumerate(self.symbols):
                    for p in TersoffPotential.twobody_parameter_names:
                         name1 = '{}{}_{}'.format(s1,s2,p)
                         name2 = '{}{}_{}'.format(s2,s1,p)
                         if name1 not in self.parameters or name2 not in self.parameters:
                             self.parameters[name1] = parameters['{}{}{}_{}'.format(s1,s2,s3,p)]
                    for p in TersoffPotential.threebody_parameter_names:
                        name1 = '{}{}{}_{}'.format(s1,s2,s3,p)
                        self.parameters[name1] = parameters[name1]

