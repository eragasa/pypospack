import copy
import yaml
from collections import OrderedDict
from pypospack.io.filesystem import OrderedDictYAMLLoader  

class PyposmatConfigurationFile(object):

    def __init__(self,filename=None):
        assert any([
            isinstance(filename,str),
            type(filename) is type(None)
            ])

        self.filename_in = None
        self.filename_out = None
        self.configuration = None

        self._parameter_names = None
        self._qoi_names = None
        self._error_names = None

        if filename is not None:
            self.read(filename=filename)

    @property
    def n_iterations(self):
        return self.sampling_type['n_iterations']
    @property
    def qois(self):
        return self.configuration['qois']

    @property
    def qoi_targets(self):
        return OrderedDict([(k,v['target']) for k,v in self.qois.items()])

    @qois.setter
    def qois(self,qois):
        assert isinstance(qois,OrderedDict)
        if self.configuration is None: self.configuration = OrderedDict()
        self.configuration['qois'] = OrderedDict()
        self.configuration['qois'] = copy.deepcopy(qois)

    @property
    def qoi_constraints(self):
        return self.configuration['qoi_constraints']
    
    @qoi_constraints.setter
    def qoi_constraints(self,qoi_constraints):
        assert isinstance(qoi_constraints,OrderedDict)
        if self.configuration is None: self.configuration = OrdereDDict()
        self.configuration['qoi_constraints'] = OrderedDict()
        self.configuration['qoi_constraints'] = copy.deepcopy(qoi_constraints)
    
    @property
    def structures(self):
        return self.configuration['structures']

    @structures.setter
    def structures(self,structures):
        assert isinstance(structures,OrderedDict)
        if self.configuration is None: self.configuration = OrderedDict()
        self.configuration['structures'] = OrderedDict()
        self.configuration['structures'] = copy.deepcopy(structures)
    
    @property
    def potential(self):
        return self.configuration['potential']

    @potential.setter
    def potential(self,potential):
        assert isinstance(potential,OrderedDict)
        if self.configuration is None: self.configuration = OrderedDict()
        self.configuration['potential'] = OrderedDict()
        self.configuration['potential'] = copy.deepcopy(potential)

    @property
    def sampling_type(self):
        return self.configuration['sampling_type']

    @sampling_type.setter
    def sampling_type(self,sampling_type):
        SAMPLING_TYPES = ['parametric','kde']
        self.configuration['sampling_type'] = copy.deepcopy(sampling_type)

    @property
    def sampling_distribution(self):
        return self.configuration['sampling_dist']

    @sampling_distribution.setter
    def sampling_distribution(self,distribution):
        self.configuration['sampling_dist'] = copy.deepcopy(distribution)

    @property
    def sampling_constraints(self):
        if 'sampling_constraints' not in self.configuration:
            self.configuration['sampling_constraints'] = None
        return self.configuration['sampling_constraints']

    @sampling_constraints.setter
    def sampling_constraints(self,sampling_constraints):
        self.configuration['sampling_constraints'] = copy.deepcopy(sampling_constraints)
        
    @property
    def mc_seed(self):
        return self.configuration['mc_seed']

    @mc_seed.setter
    def mc_seed(self,seed):
        assert type(seed) is int
        self.configuration['mc_seed'] = seed

    @property
    def parameter_distribution_definitions(self):
        return self.configuration['param_dist_def']

    @parameter_distribution_definitions.setter
    def parameter_distribution_definitions(self,param_def):
        assert isinstance(param_def,OrderedDict)
        self.configuration['param_dist_def'] = OrderedDict()
        self.configuration['param_dist_def'] = copy.deepcopy(param_def)

    @property
    def parameter_constraints(self):
        return self.configuration['param_constraints']

    @parameter_constraints.setter
    def parameter_constraints(self,constraints):
        assert isinstance(constraints,OrderedDict)
        self.configuration['param_constraints'] = OrderedDict()
        self.configuration['param_constraints'] = copy.deepcopy(constraints)

    def read(self,filename):
        self.filename_in = filename
        self.configuration = None
        with open(filename,'r') as f:
            self.configuration = yaml.load(f, OrderedDictYAMLLoader)

        self.parameter_names = [k for k in self.configuration['sampling_dist']]
        self.qoi_names = [q for q in self.configuration['qois']]
        self.error_names = ["{}.err".format(q) for q in self.qoi_names]
    
    def write(self,filename):
        self.filename_out = filename
        _configuration = copy.deepcopy(self.configuration)
        with open(filename,'w') as f:
            yaml.dump(_configuration,f, default_flow_style=False)    
