from pypospack.potential import Potential


class TersoffPotential(Potential):

    def __init__(self, symbols):
        super().__init__(symbols=symbols)
        self._pot_type = 'tersoff'
        self._determine_parameter_names()
        self._fname_potential_file = 'potential.mod'

    def _determine_parameter_names(self):
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

    def write_lammps_potential_file(self):
        for i, el1 in enumerate(self.symbols):
            for j, el2 in enumerate(self.symbols):
                for k, el3 in enumerate(self.symbols):
                    s = '{}{}{}'.format(el1, el2, el3)
                    s += '{} {} {}'.format(el1, el2, el3)
                    s += ' ' + self.parameters['{}_gamma'.format(s)]
                    s += ' ' + self.parameters['{}_lambda3'.format(s)]
                    s += ' ' + self.parameters['{}_c'.format(s)]
                    s += ' ' + self.parameters['{}_d'.format(s)]
                    s += ' ' + self.parameters['{}_costheta0'.format(s)]
                    s += ' ' + self.parameters['{}_n'.format(s)]
                    s += ' ' + self.parameters['{}_beta'.format(s)]
                    s += ' ' + self.parameters['{}_lambda2'.format(s)]
                    s += ' ' + self.parameters['{}_B'.format(s)]
                    s += ' ' + self.parameters['{}_R'.format(s)]
                    s += ' ' + self.parameters['{}_D'.format(s)]
                    s += ' ' + self.parameters['{}_lambda1'.format(s)]
                    s += ' ' + self.parameters['{}_A'.format(s)]
                    s += '\n'


