import copy
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatConfigurationFile2(PyposmatConfigurationFile):

    @property
    def structure_directory(self):
        self._structure_directory = self.structures['structure_directory']
        return self._structure_directory
    @property
    def n_iterations(self):
        self._n_iterations = self.sampling_type['n_iterations']
        return self._n_iterations

    @property
    def qoi_names(self):
        self._qoi_names = [k for k in self.qois]

    @property
    def parameter_names(self):
        self._parameter_names = [p for p in self.sampling_distribution]
        return self._parameter_names

    @property
    def parameter_constraints(self):
        if self.sampling_constraints is not None:
            self._parameter_constraints = copy.deepcopy(
                self.sampling_constriants
            )
        else:
            self._parameter_constraints = OrderedDict()

    @property
    def free_parameter_names(self):
        self._free_parameter_names = []
        for k,v in self.sampling_distribution.items():
            if v[0] != 'equals':
                self._free_parameter_names.append(k)
        return self._free_parameter_names

    @property
    def constrained_parameter_names(self):
        self._constrained_parameter_names = None
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                self._constrained_parameter_names.append(p)
        return self._constrained_parameter_names

    def print_structure_database(self):
        print(80*'-')
        print('{:^80}'.format('STRUCTURE DATABASE'))
        print(80*'-')
        print('structure_directory:{}'.format(self.structure_directory))
        print('')
        print('{:^20} {:^20}'.format('name','filename'))
        print('{} {}'.format(20*'-',20*'-'))
        for k,v in self.structures['structures'].items():
            print('{:20} {:20}'.format(k,v))

    def print_sampling_configuration(self):
        print(80*'-')
        print('{:^80}'.format('SAMPLING CONFIGURATION'))
        print(80*'-')
        print('{:^10} {:^10} {:^20}'.format(
            'iteration',
            'n_samples',
            'sampling_type'))
        print('{} {} {}'.format(10*'-',10*'-',20*'-'))

        for i in range(self.n_iterations):
            _n_samples = self.sampling_type[i]['n_samples']
            _sample_type = self.sampling_type[i]['type']
            print('{:^10} {:^10} {:^20}'.format(i,_n_samples,_sample_type))

    def print_sampling_distribution(self):
        print(80*'-')
        print('{:80}'.format('INITIAL PARAMETER DISTRIBUTION'))
        print(80*'-')
        for p in self.sampling_distribution:
            if p in self.free_parameter_names:
                str_free = 'free'
                print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
                    p,
                    str_free,
                    self.sampling_distribution[p][0],
                    self.sampling_distribution[p][1]['a'],
                    self.sampling_distribution[p][1]['b']))
            else:
                str_free = 'not_free'
                print('{:^20} {:^10}'.format(p,str_free))

if __name__ == "__main__":

    import Ni__eam__morse_exp_fs_0 as config
    configuration = PyposmatConfigurationFile2()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    configuration.write(filename='pyposmat.config.in')
    configuration.read(filename='pyposmat.config.in')

    configuration.print_sampling_configuration()
    configuration.print_sampling_distribution()
    configuration.print_structure_database()
