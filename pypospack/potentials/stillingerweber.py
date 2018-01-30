# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from pypospack.potential import Potential

class StillingerWeber(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'stillingerweber'
        self._determine_parameter_names()
        self._fname_potential_file = 'potential.mod'

    @property
    def parameter_names(self):
        return self._param_names

    def _determine_parameter_names(self):
        # TODO: This is only written for a single element potential
        symbols = self._symbols
        for i in range(n_symbols):
            for j in range(n_symbols):
                for k in range(n_symbols):
                    el1 = symbols[i]
                    el2 = symbols[j]
                    el3 = symbols[k]
                    self._add_parameter_names(el1,el2,el3)

    def _self_add_parameter_names(self,el1,el2,el3):
        s = "{}{}{}".format(el1,el2,el3)
        tersoff_param_names = ['m','gamma','lambda3','c','d','costheta0','n','beta',
                       'lambda2','B','R','D','lambda1','A']
        for p in param_names:
            self._param_names.append("{}_{}".format(s,p))

    def write_lammps_potential_file(self):
        fname_potential_mod = 'potential.mod'
        fname_tersoff_file = 'potential.tersoff'

        for i, el1 in enumerate(elements):
            for j, el2 in enumerate(elements):
                for k, el3 in enumerate(elements):
                    s = '{}{}{}'.format(el1,el2,el3)
                    m = self._param_dict['{}_m'.format(s)]
                    str_pot += '{} {} {}'.format(el1,el2,el3)
                    str_out += ' ' + self._param_dict['{}_gamma'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda3'.format(s)]
                    str_out += ' ' + self._param_dict['{}_c'.format(s)]
                    str_out += ' ' + self._param_dict['{}_d'.format(s)]
                    str_out += ' ' + self._param_dict['{}_costheta0'.format(s)]
                    str_out += ' ' + self._param_dict['{}_n'.format(s)]
                    str_out += ' ' + self._param_dict['{}_beta'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda2'.format(s)]
                    str_out += ' ' + self._param_dict['{}_B'.format(s)]
                    str_out += ' ' + self._param_dict['{}_R'.format(s)]
                    str_out += ' ' + self._param_dict['{}_D'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda1'.format(s)]
                    str_out += ' ' + self._param_dict['{}_A'.format(s)]
                    str_out += '\n'
