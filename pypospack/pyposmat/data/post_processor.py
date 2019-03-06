import os, copy
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

class PyposmatPostProcessor(object):

    def __init__(self,configuration_fn=None,datafile_fn=None):
        self.ERROR_STR_FORMAT = "{}.err"
        self.ABSOLUTE_ERROR_STR_FORMAT = "{}.abserr"
        self.NORMALIZED_ERROR_STR_FORMAT = "{}.nerr"

        self._datafile = None
        self._configuration = None
   
        self.configuration_fn = None
        self.datafile_fn = None
       
        self.qoi_fitting_names = None
        self.qoi_testing_names = None
        self.qoi_names = None

        self.error_fitting_names = None
        self.error_testing_names = None
        self.error_names = None

        if configuration_fn is not None: 
            self.read_configuration(filename=configuration_fn)
        if datafile_fn is not None: 
            self.read_datafile(filename=datafile_fn)

    @property
    def configuration(self):return self._configuration
    
    @configuration.setter
    def configuration(self,config): self._configuration = config

    @property
    def datafile(self):return self._datafile
    
    @datafile.setter
    def datafile(self,datafile): self._datafile=datafile
    
    @property
    def parameter_names(self):
        return self._parameter_names

    def read_configuration(self,filename=None):
        if filename is not None:
            self.configuration_fn = filename
        _filename = self.configuration_fn

        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)
    
        self._parameter_names = self.configuration.parameter_names

        # determine_qoi_names_from_configuration_file
        self.qoi_fitting_names = self.configuration.qoi_names
        self.qoi_testing_names = self.configuration.qoi_validation_names
        self.qoi_names = self.qoi_fitting_names + self.qoi_testing_names
        
        # determine_error_names_from_configuration_file
        self.error_fitting_names = self.configuration.error_names
        self.error_testing_names = self.configuration.qoi_validation_names
        self.error_names = self.error_fitting_names + self.error_testing_names

        # determine_absolute_error_names_from_configuration_file
        self.abs_error_fitting_names = ['{}.abserr'.format(v) for v in self.qoi_fitting_names]
        self.abs_error_testing_names = ['{}.abserr'.format(v) for v in self.qoi_testing_names]
        self.abs_error_names = self.abs_error_fitting_names + self.abs_error_testing_names

        # determine_normalized_error_names_from_configuration_file
        self.norm_error_fitting_names = ['{}.nerr'.format(v) for v in self.qoi_fitting_names]
        self.norm_error_testing_names = ['{}.nerr'.format(v) for v in self.qoi_testing_names]
        self.norm_error_names = self.norm_error_fitting_names + self.norm_error_testing_names

        if self.qoi_fitting_names is not None:
            _config_qois_fitting = self.configuration.qois
            self.qoi_fitting_targets = OrderedDict(
                [(q,_config_qois_fitting[q]['target']) for q in self.qoi_fitting_names]
            )

        if self.qoi_testing_names is not None:
            _config_qoi_testing = self.configuration.qois_validation
            self.qoi_testing_targets = OrderedDict(
                [(q,_config_qoi_testing[q]['target']) for q in self.qoi_testing_names]
            )

        self.qoi_targets = OrderedDict()
        for k,v in self.qoi_fitting_targets.items(): self.qoi_targets[k] = v
        for k,v in self.qoi_testing_targets.items(): self.qoi_targets[k] = v

        # old naming convention [DEPRECATED]
        self.qoi_validation_names = self.configuration.qoi_validation_names
        self.error_validation_names = self.configuration.qoi_validation_names

    def read_datafile(self,filename=None):
        if filename is not None:
            self.datafile_fn = filename
        _filename = self.datafile_fn
        self.datafile = PyposmatDataFile()

        self.datafile.read(filename=_filename)


        self._parameter_names = self.datafile.parameter_names
        self._qoi_names = self.datafile.qoi_names
        self._error_names = self.datafile.error_names

        self.df = copy.deepcopy(self.datafile.df)
        self.create_absolute_errors()

    def create_absolute_errors(self,qoi_names=None):
        if qoi_names is None or qoi_names is 'all':
            _qoi_names = self.qoi_names
        elif qoi_names == 'fitting':
            _qoi_names = self.qoi_fitting_names
        elif qoi_names == 'testing':
            _qoi_names = self.qoi_testing_names
        else:
            _qoi_names = list(qoi_names)
        
        for q in _qoi_names:
            aen = self.NORMALIZED_ERROR_STR_FORMAT.format(q)
            en = self.ERROR_STR_FORMAT.format(q)
            self.df[aen] = self.df[en].abs()

    def create_normalized_errors(self,normalization_type='by_qoi',qoi_names=None):
        """ normalize errors

        This class normalizes errors
        
        Args:
            normalization_type(str): normalization type
            qoi_names (str,list): allows both string arguments as well as a list of string.
        If a list of strings is provided, it will create normalized errors based on the 
        list of string.  If qoi_name is set to 'all', it will calculate normalized error based on
        both the attributes based on both the fitting qoi as well as the testing qoi.  If no
        argument is provided, normalized error will be set to 'all'.   Valid options: all, fitting, testing, or
        list of qoi strings.
            
        """

        if qoi_names is None or qoi_names is 'all':
            _qoi_names = self.qoi_names
        elif qoi_names == 'fitting':
            _qoi_names = self.qoi_fitting_names
        elif qoi_names == 'testing':
            _qoi_names = self.qoi_testing_names
        else:
            _qoi_names = list(qoi_names)

        _normalization_functions = OrderedDict([
            ('by_qoi',self.normalize_error_by_qoi)
            ])
        
        for q in _qoi_names:
            self.normalization_functions[normalization_type](qoi_names)

    def normalize_error_by_qoi(self,qoi_name,df=None):
        if df is None:
            _df = copy.deepcopy(self.df)
        elif isinstance(df,pd.DataFrame):
            nerr_name = self.NORMALIZED_ERROR_STR_FORMAT.format(q)
            err_name = self.ERROR_STR_FORMAT.format(q)
            self.df[nerr_name] = self.df[err_name]/self.qoi_targets[q] - 1.

