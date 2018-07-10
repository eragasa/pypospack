from pypospack.potential import Potential
from collections import OrderedDict
import os


class TersoffPotential(Potential):

    def __init__(self, symbols):
        super().__init__(symbols=symbols)
        self._pot_type = 'tersoff'
        self._init_parameter_names()
        self._fname_potential_file = 'potential.mod'
        self.lmps_parameter_filename = 'lmps_parameter_filename'

    def _init_parameter_names(self):
        # TODO: This is only written for a single element potential
        for i in range(len(self.symbols)):
            for j in range(len(self.symbols)):
                for k in range(len(self.symbols)):
                    el1 = self.symbols[i]
                    el2 = self.symbols[j]
                    el3 = self.symbols[k]
                    self._add_parameter_names(el1, el2, el3)

    def _add_parameter_names(self, el1, el2, el3):
        s = "{}{}{}".format(el1, el2, el3)
        tersoff_param_names = ['m','gamma','lambda3','c','d','costheta0','n','beta',
                               'lambda2','B','R','D','lambda1','A']
        for p in tersoff_param_names:
            self.param_names.append("{}_{}".format(s, p))

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

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
                    s = '{}{}{}'.format(el1, el2, el3)
                    str_out += '{} {} {}'.format(el1, el2, el3)
                    str_out += ' ' + self.parameters['{}_gamma'.format(s)]
                    str_out += ' ' + self.parameters['{}_lambda3'.format(s)]
                    str_out += ' ' + self.parameters['{}_c'.format(s)]
                    str_out += ' ' + self.parameters['{}_d'.format(s)]
                    str_out += ' ' + self.parameters['{}_costheta0'.format(s)]
                    str_out += ' ' + self.parameters['{}_n'.format(s)]
                    str_out += ' ' + self.parameters['{}_beta'.format(s)]
                    str_out += ' ' + self.parameters['{}_lambda2'.format(s)]
                    str_out += ' ' + self.parameters['{}_B'.format(s)]
                    str_out += ' ' + self.parameters['{}_R'.format(s)]
                    str_out += ' ' + self.parameters['{}_D'.format(s)]
                    str_out += ' ' + self.parameters['{}_lambda1'.format(s)]
                    str_out += ' ' + self.parameters['{}_A'.format(s)]
                    str_out += '\n'
        return str_out

