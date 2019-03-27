import os, copy, yaml
from collections import OrderedDict

import pypospack.potential
from pypospack.io.filesystem import OrderedDictYAMLLoader  

class PyposmatConfigurationFile(object):
    """Class for reading, writing, and analyzing configuration file

    Args:
        filename(str,optional): the filename of the configuration file.  This
           is the location which all IO will take place with.  If this
           argument is passed, then the class will automatically read
           the class

    Notes:
        To read a file:
        
        from pypospack.pyposmat import PyposmatConfigurationFile
        o = PyposmatConfigurationFile(filename=_filename)

        or

        from pypospack.pyposmat import PyposmatConfigurationfile
        o = PyposmatConfigurationFile()
        o.read(filename=_filename)


    """

    def __init__(self,filename=None):

        # type checking the arguments
        assert type(filename) in [type(None),str]
      
        # public attributes
        self.filename = None #: str: the name for IO

        self.filename_in = None
        self.filename_out = None
        self.configuration = None

        # private attributes
        self._parameter_names = None
        self._qoi_names = None
        self._error_names = None

        if filename is str:
            if os.path.isabs(filename): 
                self.filename = filename
            else:
                self.filename = os.path.abspath(filename)
            
        if self.filename is str:
            self.read(filename=self.filename)

    @property
    def n_iterations(self):
        """int:the number of iterations in this fitting process."""
        return self.sampling_type['n_iterations']

    @property
    def qois(self):
        return self.configuration['qois']

    @qois.setter
    def qois(self,qois):
        assert isinstance(qois,OrderedDict)
        if self.configuration is None: self.configuration = OrderedDict()
        self.configuration['qois'] = OrderedDict()
        self.configuration['qois'] = copy.deepcopy(qois)
   
    @property
    def qoi_fitting(self):
        return self.configuration['qoi_f']
    @property
    def qoi_fitting_values(self):
        return OrderedDict([(k,v['target']) for k,v in self.qoi_fitting.items()])

    @property
    def qoi_testing(self):
        return self.configuration['qoi_t']

    @property
    def qoi_testing_values(self):
        return OrderedDict([(k,v['target']) for k,v in self.qoi_testing.items()])

    @property
    def qoi_targets(self):
        return OrderedDict([(k,v['target']) for k,v in self.qois.items()])
   
    @property
    def qoi_validation_targets(self):
        try:
            return OrderedDict([(k,v['target']) for k,v in self.qois_validation.items()])
        except KeyError as e:
            return None

    @property
    def qois_validation(self):
        try:
            return self.configuration['qois_v']
        except KeyError as e:
            return None

    @qois_validation.setter
    def qois_validation(self,qois):
        assert isinstance(qois,OrderedDict)
        if self.configuration is None:
            self.configuration = OrderedDict()

        self.configuration['qois_v'] = OrderedDict()
        self.configuration['qois_v'] = copy.deepcopy(qois)

    @property
    def qoi_validation_targets(self):
        if type(self.qois_validation) is type(None):
            return None
        else:
            return OrderedDict([(k,v['target']) for k,v in self.qois_validation.items()])

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

    # this property is broken
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

    @property
    def reference_potentials(self):
        if 'reference_potentials' not in self.configuration:
            return None
        else:
            return self.configuration['reference_potentials']

    @reference_potentials.setter
    def reference_potentials(self,potentials):
        self.configuration['reference_potentials'] = copy.deepcopy(potentials)

    @property
    def latex_labels(self):
        try:
            return self.configuration['latex_labels']
        except KeyError as e:
            if e.args[0] == 'latex_labels':
               self.configuration['latex_labels'] = OrderedDict()
               for v in self.qoi_names:
                   self.configuration['latex_labels'][v] = v
               for v in self.parameter_names:
                   self.configuration['latex_labels'][v] = v
               for v in self.error_names:
                   self.configuration['latex_labels'][v] = v

               return self.configuration['latex_labels']

    @latex_labels.setter
    def latex_labels(self,labels):
        assert isinstance(labels,OrderedDict)
        self.configuration['latex_labels'] = OrderedDict()
        self.configuration['latex_labels'] = copy.deepcopy(labels)

    def read(self,filename):
        self.filename_in = filename
        self.configuration = None
        with open(filename,'r') as f:
            self.configuration = yaml.load(f, OrderedDictYAMLLoader)

        self.set_parameter_names_from_configuration_dictionary()
        self.set_free_parameter_names_from_configuration_dictionary()
        
        # set qoi_names,error_names, and normalized_error_names
        _qois = self.configuration['qois']
        self.qoi_names = [q for q in _qois]
        self.error_names = ["{}.err".format(q) for q in _qois]
        self.normed_error_names = ["{}.nerr".format(q) for q in _qois]

        try:
            # set qoi_validation,error_validation_names, and normalized_error_names if
            # the field exists
            _qoi_v_names = self.configuration['qois_v']
            self.qoi_validation_names = [q for q in _qoi_v_names]
            self.error_validation_names = ["{}.err".format(q) for q in _qoi_v_names]
            self.normed_error_validation_names = ['{}.nerr'.format(q) for q in _qoi_v_names]
        except KeyError as e:
            # set up attribute to default NoneType if they aren't specified
            self.qoi_validation_names = None
            self.error_validation_names = None
            self.normed_error_validation_names = None

    def set_parameter_names_from_configuration_dictionary(self,o_config=None):

        # o_config is an OrderedDictionary object.
        # o_config is read in from a YAML file and stored as the self.config
        if o_config is None:
            _o_config = self.configuration

        # the keys of the 'sampling_dist' are the parameter_names, they need to match
        # the naming convention defined in the Potential
        self.parameter_names = [k for k,v in _o_config['sampling_dist'].items()]

        return self.parameter_names

    def set_free_parameter_names_from_configuration_dictionary(self,o_config=None):
        """
        Args:
            o_config (collections.OrderedDict): a configuration dictionary which was
                read in using YAML for marshalling and unmarshalling of the object.  By
                default this is set to None, which uses the existing configuration dictionary
        """
        # o_config is an OrderedDictionary object.
        # o_config is read in from a YAML file and stored as the self.config
        if o_config is None:
            _o_config = self.configuration

        # the key object is the name of the parameter
        # the value object is a list
        # v[0]:
        #   'equals' = equality constraint on the parameter defined in v[1]
        #   'uniform' = uniform distribution with lower bounds and upper bounds defined in v[1]
        #       dictionary object with v[1]['a'] being the lower bound, and v[1]['b'] being the
        #       upper bound
        #   'normal' = normal distribution with the mean and standard deviation defined in 
        #       v[1] with v[1]['mu'] being the mean and v[1]['sigma'] being the standard deviation.
        self.free_parameter_names = \
                [k for k,v in _o_config['sampling_dist'].items() if v[0] != 'equals']

        return self.parameter_names

    def write(self,filename):
        self.filename_out = filename
        _configuration = copy.deepcopy(self.configuration)
        with open(filename,'w') as f:
            yaml.dump(_configuration,f, default_flow_style=False)

    def str__structure_database(self):
        _structure_dir = self.structures['structure_directory']
        _structures = self.structures['structures']
        str_out = [
            80*'-',
            '{:^80}'.format('STRUCTURE DATABASE'),
            80*'-',
            'structure_directory:{}'.format(_structure_dir),
            '',
            '{:^20} {:^20}'.format('name','filename'),
            '{} {}'.format(20*'-',20*'-')
        ]
        for k,v in _structures.items():
            str_out.append('{:20} {:20}'.format(k,v))

        return "\n".join(str_out)

    # figures for
    def get_ok_fail(self,b):
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'

        if b:
            s = '['+OKGREEN+'OK'+ENDC+']  '
        else:
            s = '['+FAIL+'FAIL'+ENDC+']'

        return s

    def validate(self):
        self.validate_potential(potential=self.potential)
        self.validate_parameters(potential=self.potential,parameters=self.parameter_names)
        self.validate_structure_db(structures=self.structures)
        self.validate_sampling_type(sampling_type=self.sampling_type)

    def validate_sampling_type(self,sampling_type):
        _n_iterations = sampling_type['n_iterations']
        if type(_n_iterations) is not int:
            s = ['n_iterations must be int, {}'.format(_n_iterations)]
            print("\n".join(s))
        else:
            s = "{:5} n_iterations:{}".format(self.get_ok_fail(True),_n_iterations)
            print(s)

        _mc_seed = sampling_type['mc_seed']
        if _mc_seed is None:
            s = '{:5} mc_seed:{}, using automatic seed generation'.format(
                    self.get_ok_fail(True),
                    _mc_seed
                )
            print(s)
        elif type(_mc_seed) is int:
            s = '{:5} mc_seed:{}'.format(
                    self.get_ok_fail(True),
                    _mc_seed
                )
            print(s)

        if all([v in sampling_type for v in range(_n_iterations)]):
            s = [
                "{:5} n_iterations of sampling_types defined".format(
                    self.get_ok_fail(True))
                ]
            print("\n".join(s))
        else:
            s = [
                "{:5} not all iterations have been defined".format(
                    self.get_ok_fail(False))
                ]
            print("\n".join(s))
        
        for i in range(_n_iterations):
            _sampling_type = sampling_type[i]['type']
            if _sampling_type == 'parametric':
                _is_supported_sampling_type = True

                try:
                    _n_samples = sampling_type[i]['n_samples']
                except NameError as e:
                    s = "{:5} n_samples must be defined for {} sampling_type in iteration {}".format(
                        self.get_ok_fail(False),
                        _sampling_type,
                        i)
                    print(s)
                    _is_supported_sampling_type = False

                if _is_supported_sampling_type:
                    s = "{:5} iteration {}: {} sampling".format(
                        self.get_ok_fail(True),
                        i,
                        _sampling_type
                    )
                    print(s)

            elif _sampling_type == 'kde':
                _is_supported_sampling_type = True
                try:
                    _n_samples = sampling_type[i]['n_samples']
                except NameError as e:
                    s = "{:5} n_samples must be defined for {} sampling_type in iteration {}".format(
                        self.get_ok_fail(False),
                        i,
                        _sampling_type
                    )
                    print(s)
                if _is_supported_sampling_type:
                    s = "{:5} iteration {}: {} sampling".format(
                        self.get_ok_fail(True),
                        i,
                        _sampling_type
                    )
                    print(s)

            elif _sampling_type == 'from_file':
                _is_supported_sampling_type = True
                if _is_supported_sampling_type:
                    s = "{:5} iteration {}: {} sampling".format(
                        self.get_ok_fail(True),
                        i,
                        _sampling_type)
                    print(s)
            else:
                _is_supported_sampling_type = False
    
    def validate_structure_db(self,structures):
        _structure_dir = self.structures['structure_directory']
        _structure_dir_exists = os.path.isdir(_structure_dir)
        _structures = self.structures['structures']
        
        if _structure_dir_exists:
            s = ['structure_directory:{}'.format(_structure_dir)]
        else:
            s = ['structure_directory does not exist, {}'.format(_structure_dir)]
        
        _all_structures_exist = True
        for _structure_name,_structure_fn in _structures.items():
            _structure_fn = os.path.join(_structure_dir,_structure_fn)                        
            _structure_exists = os.path.isfile(_structure_fn)

            if not _structure_exists:
                _all_structures_exist = False
            s += ['{:10} {:20} {}'.format(self.get_ok_fail(_structure_exists),_structure_name,_structure_fn)]

        print("\n".join(s))

    def validate_potential(self,potential):
        if 'potential_type' not in potential:
            s = "configuration file requires a potential type"
            print(s)
            _is_supported_potential=False
        else:
            _potential_type = potential['potential_type']
            
            _is_supported_potential = self.is_supported_potential(potential_type=_potential_type)

            if _is_supported_potential:
                s = "potential_type:{}".format(_potential_type)
                print(s)
                if _potential_type == 'eam':
                    _is_supported_potential = self.validate_eam_potential(potential)
                else:
                    _is_supported_potential = True
            else:
                s = "potential_type not supported:{}".format(_potential_type)
                print(s)
                _is_supported_potential = True


    def is_supported_potential(self,potential_type):
        return potential_type in pypospack.potential.get_supported_potentials()

    def validate_eam_potential(self,potential):
        _func_pair = potential['pair_type'] 
        _func_density = potential['density_type']
        _func_embedding = potential['embedding_type']

        _is_supported_potential = True
        if not _func_pair in pypospack.potential.pair_potentials:
            s = "pair_type is not supported:{}".format(_func_pair)
            print(s)
            _is_supported_potential = False
        else:
            s = "pair_type:{}".format(_func_pair)
            print(s)

        if not _func_density in pypospack.potential.eam_density_functions:
            s = "'density_type is not supported:{}".format(_func_density)
            print(s)
            _is_supported_potential = False
        else:
            s = "density_type:{}".format(_func_density)

        if not _func_embedding in pypospack.potential.eam_embedding_functions:
            s = "embedding_type is not supported:{}".format(_func_embedding)
            print(s)
            _is_supported_potential = False
        else:
            s = "embedding_type:{}".format(_func_embedding)
            print(s)

        return _is_supported_potential

    def validate_parameters(self,potential,parameters):
        _potential_type = potential['potential_type']
        _symbols = potential['symbols']

        if _potential_type == 'eam':
            # EAM potentials have a different constructor interface different than
            # other potentials
            _pair_type = potential['pair_type']
            _density_type = potential['density_type']
            _embedding_type = potential['embedding_type']

            _obj_potential = pypospack.potential.EamPotential(
                    symbols=_symbols,
                    func_pair=_pair_type,
                    func_density=_density_type,
                    func_embedding=_embedding_type)
        else:
            # using the default constructor using dynamic module and class
            # loading for all other potentials
            _module_name, _class_name = pypospack.potential.PotentialObjectMap(potential_type='all')
            _module = importlib.import_module(_module_name)
            _class = getattr(_module,_class_name)
            _obj_potential = _class(symbols=_symbols)

        # checking if all parameters for the potential are defined
        s = ['checking if required parameters are defined']
        all_parameters = set.union(
                set(parameters),
                set(_obj_potential.parameter_names)
        )
        for pn in all_parameters:
            pn_in_parameter_set = pn in parameters
            pn_in_potential_definition = pn in _obj_potential.parameter_names
            pn_is_ok = pn_in_parameter_set and pn_in_potential_definition
            s.append(
                '{:8} {:20} {:6} {:6}'.format(
                    self.get_ok_fail(pn_is_ok),
                    pn,
                    self.get_ok_fail(pn_in_parameter_set),
                    self.get_ok_fail(pn_in_potential_definition)
                )
            )

        _all_parameters_validated = True
        s.append('{:5} {}'.format(
            self.get_ok_fail(set(parameters)==set(_obj_potential.parameter_names)),
            'parameters passed are the same as the parameters required'
        ))
        for _parameter_name in _obj_potential.parameter_names:
            s.append('{:5} {:15}'.format(
                self.get_ok_fail(_parameter_name in parameters),
                _parameter_name)
            )
            if _parameter_name not in parameters:
                _all_parameters_validated = False
       
        s.append(
                'parameters in parameter set'
        )

        for pn in parameters:
            s.append(
                '{:5} {:15}'.format(
                    self.get_ok_fail(pn in _obj_potential.parameter_names),
                    pn
                )
            )

        print('\n'.join(s))
        return _all_parameters_validated

    def get_free_parameter_names(self):
        for pn in self.parameter_names:
            print(pn)

class PyposmatConfigurationFileValidator(object):
    pass
