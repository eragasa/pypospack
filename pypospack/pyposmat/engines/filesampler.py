import os
import pandas as pd
import numpy as np
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.engines import PyposmatEngine

class PyposmatFileSampler(PyposmatEngine):

    def __init__(self,
            config_fn='pyposmat.config.in',
            data_in_fn='pyposmat.results.in',
            data_out_fn='pyposmat.results.out',
            mpi_rank = None,
            mpi_size=None,
            log = None,
            base_directory = None):
        
        PyposmatEngine.__init__(self,
                filename_in=config_fn,
                filename_out=data_out_fn,
                base_directory=base_directory,
                fullauto=False)

        self.configuration_fn = config_fn
        self.datafile_in_fn = data_in_fn
        self.datafile_out_fn = data_out_fn

        self.mpi_rank = mpi_rank
        self.mpi_size = mpi_size

        self.configuration = None
        self.datafile = None
        self.subselect_df = None
        self.reference_potentials = None

        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None
        self.normed_error_names = None
        
        self.qoi_validation_names = None
        self.error_validation_names = None
        self.normed_error_validation_names = None

        self.qoi_targets = None
        self.qoi_validation_target = None

        if config_fn is not None:
            self.read_configuration_file(config_fn=config_fn)

        if data_in_fn is not None:
            self.read_datafile_in(datafile_fn=data_in_fn)

        if data_out_fn is not None:
            self.datafile_out = PyposmatDataFile(data_out_fn)

    def __log(self,m):
        print(m)

    def configure_qoi_manager(self,qois=None,use_fitting_qois=True,use_testing_qois=False):

        if qois is not None:
            _qois = copy.deepcopy(qois)
        else:
            _qois = OrderedDict()

            if use_fitting_qois:
                for k,v in self.configuration.qois.items():
                    _qois[k]=v

            if use_testing_qois:
                for k,v in self.configuration.qois_validation.items():
                    _qois[k]=v
        PyposmatEngine.configure_qoi_manager(self,_qois)

    def configure_task_manager(self):
        PyposmatEngine.configure_task_manager(self)

    def read_configuration_file(self,config_fn=None):
        if config_fn is not None:
            self.configuration_fn = config_fn
        
        _config_fn = self.configuration_fn


        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=_config_fn)

        self.structure_directory = self.configuration.structures['structure_directory']
        if os.path.isdir(self.structure_directory):
            msg = "[OK] structure_directory:".format(self.structure_directory)
            self.__log(msg)
        else:
            msg = "[FAIL] structure_directory:".format(self.structure_directory)
            raise PyposmatEngineError(msg)
        
        # set parameter names
        self.parameter_names = self.configuration.parameter_names
        
        # set name arrays for testing qois
        self.qoi_names = self.configuration.qoi_names
        self.error_names = self.configuration.error_names
        self.normed_error_names = self.configuration.normed_error_names

        # set name arrays for validation qois
        self.qoi_validation_names = self.configuration.qoi_validation_names
        self.error_validation_names = self.configuration.error_validation_names
        self.normed_error_validation_names = self.configuration.normed_error_validation_names

        # set dictionaries for qoi targets
        self.qoi_targets = self.configuration.qoi_targets
        self.qoi_validation_targets = self.configuration.qoi_validation_targets
    
        # set dictionary for reference potentials
        self.reference_potentials = self.configuration.reference_potentials

    def read_datafile_in(self,datafile_fn):
        self.datafile_in = PyposmatDataFile()
        self.datafile_in.read(filename=datafile_fn)

    def configure_datafile_out(self,datafile_fn=None):
        if datafile_fn is not None:
            self.datafile_out_fn = datafile_fn

        _datafile_fn = self.datafile_out_fn

        self.datafile_out = PyposmatDataFile(_datafile_fn)

    def subselect_by_dmetric(self,nsmallest=50):

        # calculated normalized errors for qois
        for iqn,qn in enumerate(self.qoi_names):
            error_name = "{}.err".format(qn)
            normed_error_name = "{}.nerr".format(qn)

            q = self.qoi_targets[qn]
            error = self.datafile_in.df[error_name]

            self.datafile_in.df[normed_error_name] = error/q

        self.datafile_in.df['d_metric'] = np.sqrt(np.square(
            self.datafile_in.df[self.normed_error_names]).sum(axis=1))  
    
        self.subselect_df = self.datafile_in.df.nsmallest(nsmallest,'d_metric')

        return self.subselect_df
    
    def run_simulations(self):

        if self.qoi_validation_names is not None:
            self.datafile_out.write_header_section(
                filename = self.datafile_out_fn,
                parameter_names = self.parameter_names,
                qoi_names = self.qoi_names,
                error_names = self.error_names,
                qoi_v_names = self.qoi_validation_names,
                error_v_names = self.error_validation_names
                )
            
        if self.reference_potentials is not None:
            self._sample_from_reference_potentials()

        if self.subselect_df is not None:
            self._sample_from_subselect_df(
                    subselect_df=self.subselect_df)
        else:
            self._sample_from_subselect_df(
                    subselect_df=self.datafile_in.df)

    def _sample_from_reference_potentials(self,reference_potentials=None):
        """

        This method assumes that the reference potentials have the same functional form as 
        the potentials being tested.

        """
        if reference_potentials is None:
            _rpotentials = self.reference_potentials

        for potential_name,potential in _rpotentials.items():
            
            try:
                _sim_id = int(float(potential_name))
            except ValueError as e:
                _sim_id = potential_name

            parameters = potential['parameters']
            
            evals = self.evaluate_parameter_set(parameters=parameters)

            _results = OrderedDict()
            _results['parameters'] = parameters
            _results['qois'] = OrderedDict()
            for v in self.qoi_names:
                try:
                    _results['qois'][v] = potential['qoi'][v]
                except KeyError as e:
                    if v in evals['qois']:
                        _results['qois'][v] = evals['qois'][v]
                    else:
                        _results['qois'][qn] = np.NaN
            _results['errors'] = OrderedDict()
            for v in self.error_names:
                try:
                    qn = ".".join([s for s in v.split(".") if s != 'err'])
                    qhat = potential['qoi'][qn]
                    q = self.configuration.qois[qn]['target']
                    _results['errors'][v] = qhat - q
                except KeyError as e:
                    if v in evals['errors']:
                        _results['errors'][v] = evals['errors'][v]
                    else:
                        _results['errors'][qn] = np.NaN

            _results['qois_v'] = OrderedDict()
            for v in self.qoi_validation_names:
                if v in evals['qois']:
                    _results['qois_v'][v] = evals['qois'][v]
                else:
                    _results['qois_v'][qn] = np.NaN
                    
            _results['errors_v'] = OrderedDict()
            for v in self.error_validation_names:
                if v in evals['errors']:
                    _results['errors_v'][v] = evals['errors'][v]
                else:
                    _results['errors_v'][qn] = np.NaN
          
            print('simulation_finished, sim_id:{}'.format(_sim_id))
            self.datafile_out.write_simulation_results(
                    sim_id = _sim_id,
                    results = _results)


    def _sample_from_subselect_df(self,subselect_df=None):
        
        if subselect_df is None:
            _subselect_df = self.subselect_df
        else:
            _subselect_df = subselect_df

        for idx,row in _subselect_df.iterrows():
            parameters = OrderedDict([(pn,row[pn]) for pn in self.parameter_names])
           
            try:
                _sim_id = int(float(row['sim_id']))
            except ValueError as e:
                _sim_id = potential_name

            evals = self.evaluate_parameter_set(parameters=parameters)
            _results = OrderedDict()
            _results['parameters'] = parameters
            _results['qois'] = OrderedDict()
            for v in self.qoi_names:
                try:
                    _results['qois'][v] = row[v]
                except KeyError as e:
                    if v in evals['qois']:
                        _results['qois'][v] = evals['qois'][v]
                    else:
                        _results['qois'][qn] = np.NaN
                    
            _results['errors'] = OrderedDict()
            for v in self.error_names:
                try:
                    _results['errors'][v] = row[v]
                except KeyError as e:
                    if v in evals['errors']:
                        _results['errors'][v] = evals['errors'][v]
                    else:
                        _results['errors'][qn] = np.NaN

            _results['qois_v'] = OrderedDict()
            for v in self.qoi_validation_names:
                try:
                    _results['qois_v'][v] = row[v]
                except KeyError as e:
                    if v in evals['qois']:
                        _results['qois_v'][v] = evals['qois'][v]
                    else:
                        _results['qois_v'][qn] = np.NaN
                    
            _results['errors_v'] = OrderedDict()
            for v in self.error_validation_names:
                try:
                    _results['errors_v'][v] = row[v]
                except KeyError as e:
                    if v in evals['errors']:
                        _results['errors_v'][v] = evals['errors'][v]
                    else:
                        _results['errors_v'][qn] = np.NaN
           
            print('simulation_finished, sim_id:{}'.format(_sim_id))
            self.datafile_out.write_simulation_results(
                    sim_id = _sim_id,
                    results = _results)

if __name__ == "__main__":
    pass
