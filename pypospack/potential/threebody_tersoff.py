from pypospack.potential import ThreeBodyPotential
from collections import OrderedDict
import os


class TersoffPotential(ThreeBodyPotential):

    potential_type = 'tersoff'
    threebody_parameter_names = ['m','gamma','lambda3','c','d','costheta0','n','beta',
                           'lambda2','B','R','D','lambda1','A']
    def __init__(self, symbols):
        potential_type = TersoffPotential.potential_type
        ThreeBodyPotential.__init__(self,
                                    symbols=symbols,
                                    potential_type=potential_type)
        self.lmps_parameter_filename = 'lmps_parameter_filename'

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

        str_out += "pair_style tersoff\n"

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

