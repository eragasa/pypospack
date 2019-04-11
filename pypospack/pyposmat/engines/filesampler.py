import os
import time
import pandas as pd
import numpy as np
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatLogFile
from pypospack.pyposmat.engines import PyposmatEngine

# necessary errors which this class needs to handle
from pypospack.exceptions import LammpsSimulationError
from pypospack.exceptions import PyposmatBadParameterError
from pypospack.exceptions import PypospackBadKdeBandwidthType
from pypospack.exceptions import PypospackTaskManagerError 


class PyposmatFileSampler(PyposmatEngine):
    """ samples from a datafile

    Args:
        config_fn(str): path of the configuration file
        data_in_fn(str): path of the datafile to sample from
        data_out_fn(str): path of the datafile to write results to
        mpi_rank(int): the MPI rank of this process
        mpi_rank(int): the MPI size of this process
        log(PyposmatLogFile): log object instance

    Attributes:
        configuration(PyposmatConfigurationFile): configuration file
        datafile_in(PyposmatDataFile): object instance of the datafile being read in
        datafile_out(PyposmatDataFile): object instance of the datafiel being written to

    """

    def __init__(self,
            config_fn='pyposmat.config.in',
            data_in_fn='pyposmat.results.in',
            data_out_fn='pyposmat.results.out',
            mpi_rank = None,
            mpi_size=None,
            o_log = None,
            log_to_stdout=True,
            base_directory = None,
            fullauto=True):

        self.DEBUG = False
        self.configuration = None
        self.datafile_in = None
        self.datafile_out = None

        self.mpi_rank = None
        self.mpi_size = None

        self.qoi_mananger = None
        self.task_manager = None

        PyposmatEngine.__init__(self,
                filename_in=config_fn,
                filename_out=data_out_fn,
                base_directory=base_directory,
                fullauto=fullauto)

        self.configuration_fn = config_fn
        self.datafile_in_fn = data_in_fn
        self.datafile_out_fn = data_out_fn

        self.configuration = None
        self.datafile = None
        self.subselect_df = None
        self.reference_potentials = None

        self.qoi_validation_names = None
        self.error_validation_names = None
        self.normed_error_validation_names = None

        self.qoi_validation_target = None
        self.obj_log = None
        self.log_to_stdout = None

        self.initialize_mpi_information(mpi_rank=mpi_rank,
                                        mpi_size=mpi_size)
        self.configure_logger(o_log=o_log,
                              log_to_stdout=log_to_stdout)
        if fullauto is True:
            self.read_configuration_file(filename=self.configuration_fn)
            self.read_datafile_in(filename=data_in_fn)
            self.configure_datafile_out(filename=data_out_fn)

    @property
    def n_iterations(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.sampling_type['n_iterations']
        else:
            return None

    @property
    def parameter_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.parameter_names
        else:
            return None

    @property
    def qoi_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.qoi_names
        else:
            return None

    @property
    def qoi_targets(self):
        if isntance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.qoi_targets
        else:
            return None

    @property
    def error_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.error_names
        else:
            return None

    @property
    def normed_error_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.normed_error_names
        else:
            return None
    
    @property
    def parameter_distribution_definition(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.sampling_distribution
        else:
            return None

    @property
    def free_parameter_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.free_parameter_names
        else:
            return None

    @property
    def parameter_constraints(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.sampling_constraints
        else:
            return None

    def initialize_mpi_information(self,mpi_rank=None,mpi_size=None):
        if isinstance(mpi_rank, int) and isinstance(mpi_size,int):
            self.mpi_rank = mpi_rank
            self.mpi_size = mpi_size
        elif mpi_rank is None and mpi_size is None:
            self.mpi_rank = 0
            self.mpi_size = 1
        else:
            error_message = "mpi_rank:{}\n".format(mpi_rank)
            error_message += "mpi_size:{}".format(mpi_size)
            raise TypeError(error_message)


    def log(self,str_msg):
        """ log message

        Args:
            str_msg(str,list): This is the message to log.  If this argument is a string, then
                the string will bemessage to log

        """
        if type(str_msg) is str:
            m = str_msg
        elif type(str_msg) is list:
            m = "\n".join(str_msg)

        if type(self.obj_log) is PyposmatLogFile:
            self.obj_log.write(m)

        if self.log_to_stdout:
            print(m)
    
    def configure_logger(self,o_log=None,log_to_stdout=True):
        """
        Configurtion of the log object has different behavior based upon the type passed
        into the argument o_log.  If o_log is PyposmatLogFile, that object will be accessed
        by reference.  A string is assumed to be a filename location.  By default the
        argument for o_log is None, which means logging will go to standard out by means of 
        the print() function.

        Args:
            o_log (str,PyposmatLogFile,None): default: None
        """

        if type(o_log) is PyposmatLogFile:
            self.obj_log = o_log
        elif type(o_log) is str:
            self.obj_log = PyposmatLogFile(filename=o_log)
        elif o_log is None:
            self.obj_log = None
        else:
            m = "log object must be str, PyposmatLogFile, or None"
            raise TypeError(m)

        if isinstance(log_to_stdout,bool):
            self.log_to_stdout = log_to_stdout
        else:
            m = "log_to_stdout must be boolean"
            raise TypeError(m)
    
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

    def read_configuration_file(self,filename=None):
        PyposmatEngine.read_configuration_file(self,filename=filename)

        # does this have to be removed?
        # self.structure_directory = self.configuration.structures['structure_directory']
        if self.DEBUG:
            if os.path.isdir(self.structure_directory):
                msg = "[OK] structure_directory:".format(self.structure_directory)
                self.__log(msg)
            else:
                msg = "[FAIL] structure_directory:".format(self.structure_directory)
                raise PyposmatEngineError(msg)
        

        # set name arrays for validation qois
        self.qoi_validation_names = self.configuration.qoi_validation_names
        self.error_validation_names = self.configuration.error_validation_names
        self.normed_error_validation_names = self.configuration.normed_error_validation_names

        # set dictionaries for qoi targets
        self.qoi_validation_targets = self.configuration.qoi_validation_targets
    
        # set dictionary for reference potentials
        self.reference_potentials = self.configuration.reference_potentials

    def read_datafile_in(self,filename=None):
        self.datafile_in = PyposmatDataFile()
        self.datafile_in.read(filename=filename)

    def configure_datafile_out(self,filename=None):
        if filename is not None:
            self.datafile_out_fn = filename

        self.datafile_out = PyposmatDataFile(self.datafile_out_fn)

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
    
    def run_simulations(self,i_iteration,n_samples=None,filename=None):
        """ run simulations

        Args:
            i_iterations(int): the iteration cycle we are on
            n_samples(int,optional): does not do anything
            filename(str,optional): the filename we are sampling from

        """
        if isinstance(filename,str):
            self.read_datafile_in(filename=filename)
        else:
            _filename = self.configuration.sampling_type[i_iteration]['file']
            if os.path.isabs(_filename):
                self.read_datafile_in(filename=_filename)
            else:
                self.read_datafile_in(filename=_filename)

        if self.qoi_validation_names is not None:
            self.datafile_out.write_header_section(
                filename = self.datafile_out_fn,
                parameter_names = self.parameter_names,
                qoi_names = self.qoi_names,
                error_names = self.error_names,
                qoi_v_names = self.qoi_validation_names,
                error_v_names = self.error_validation_names
            )
        else:
            self.datafile_out.write_header_section(
                filename = self.datafile_out_fn,
                parameter_names = self.parameter_names,
                qoi_names = self.qoi_names,
                error_names = self.error_names,
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
          
            self.datafile_out.write_simulation_results(
                    sim_id = _sim_id,
                    results = _results)


    def _sample_from_subselect_df(self,subselect_df=None):
        
        if subselect_df is None:
            _subselect_df = self.subselect_df
        else:
            _subselect_df = subselect_df

        time_start_iteration = time.time()
        n_errors = 0

        for i_sample,row in _subselect_df.iterrows():
            if i_sample%self.mpi_size != self.mpi_rank:
                pass
            else:

                # populate local parameter dictionary
                parameters = OrderedDict([(pn,row[pn]) for pn in self.parameter_names])
               
                try:
                    sim_id = int(float(row['sim_id']))
                except ValueError as e:
                    sim_id = row['sim_id']

                try:
                    evals = self.evaluate_parameter_set(parameters=parameters)
                except PyposmatBadParameterError as e:
                    n_errors += 1
                except LammpsSimulationError as e:
                    n_errors += 1
                except PypospackTaskManagerError as e:
                    n_errors += 1
                else:
                    results = OrderedDict()
                    results['parameters'] = parameters
                    results['qois'] = OrderedDict([(v,evals['qois'][v]) for v in self.qoi_names])
                    results['errors'] = OrderedDict([(v,evals['errors'][v]) for v in self.error_names])
                    self.datafile_out.write_simulation_results(sim_id = sim_id,results = results)
                finally:
                    if (i_sample+1)%10 == 0:
                        n_samples_completed = i_sample+1
                        time_end = time.time()
                        time_total = time_end-time_start_iteration
                        avg_time = time_total/n_samples_completed
                        s = 'R{}:{} samples completed in {:4f}s, Avg_time={:4f}. n_errors={}'
                        s = s.format(self.mpi_rank,n_samples_completed,time_total,avg_time,n_errors)
                        self.log(s)
    
    def print_structure_database(self):
        m = [
            80*'-',
            '{:^80}'.format('STRUCTURE DATABASE'),
            80*'-',
            'structure_directory:{}'.format(self.structure_directory),
            '',
            '{:^20} {:^20}'.format('name','filename'),
            '{} {}'.format(20*'-',20*'-')
        ]
        m += ['{:20} {:20}'.format(k,v) for k,v in self.structures['structures'].items()]
        self.log(m)

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
            _sample_type = self.configuration.sampling_type[i]['type']
            if _sample_type == 'kde_w_clusters':
                _n_samples = self.configuration.sampling_type[i]['n_samples_per_cluster']
            else:
                _n_samples = self.configuration.sampling_type[i]['n_samples']
            print('{:^10} {:^10} {:^20}'.format(i,_n_samples,_sample_type))
    
    def print_initial_parameter_distribution(self):
        self.log(80*'-')
        self.log('{:80}'.format('INITIAL PARAMETER DISTRIBUTION'))
        self.log(80*'-')
        for p in self.parameter_distribution_definition:
            if p in self.free_parameter_names:
                str_free = 'free'
                if self.parameter_distribution_definition[p][0] == 'uniform':
                    print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
                        p,
                        str_free,
                        self.parameter_distribution_definition[p][0],
                        self.parameter_distribution_definition[p][1]['a'],
                        self.parameter_distribution_definition[p][1]['b']))
                elif self.parameter_distribution_definition[p][0] == 'normal':
                    print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
                        p,
                        str_free,
                        self.parameter_distribution_definition[p][0],
                        self.parameter_distribution_definition[p][1]['mu'],
                        self.parameter_distribution_definition[p][1]['sigma']))
                else:
                    _distribution_type = self.parameter_distribution_defintion[p][0]
                    s = "incorrection parameter distribution for parameter {}.  probability distribution function, {}, is not supported"
                    s = s.format(p,_distribution_type)
                    raise ValueError(s)

            else:
                str_free = 'not_free'
                print('{:^20} {:^10}'.format(p,str_free))
if __name__ == "__main__":
    pass
