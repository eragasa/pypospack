from pypospack.potential import Potential

class ThreeBodyPotential(Potential):

    def _init_parameter_names(self):
        # TODO: This is only written for a single element potential
        for i in range(len(self.symbols)):
            for j in range(len(self.symbols)):
                for k in range(len(self.symbols)):
                    el1 = self.symbols[i]
                    el2 = self.symbols[j]
                    el3 = self.symbols[k]
                    self._add_parameter_names(el1, el2, el3)

    def _add_parameter_names(self,el1,el2,el3):
        s = "{}{}{}".format(el1,el2,el3)

        self.parameter_names = []
        for p in self.threebody_parameter_names:
            self.parameter_names.append('{}__{}'.format(s,p))

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None
