import os,shutil,sys
import numpy as np
from mpi4py import MPI
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
# from pypospack.pyposmat.engines import PyposmatMonteCarloSampler
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data.cluster_analysis import PyposmatClusterAnalysis
# ---- imports for PyposmatMonteCarloSampler
import time,sys,os,copy,shutil,importlib
from collections import OrderedDict
import numpy as np
import scipy.stats
from pypospack.kde import Chiu1999_h
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.task.lammps import LammpsSimulationError
from pypospack.task.task_manager import PypospackTaskManagerError
from pypospack.potential import PotentialObjectMap
from numpy.linalg import LinAlgError

# ---- additional imports for PyposmatMonteCarloSampler
from pypospack.qoi import QoiManager
from pypospack.task import TaskManager
from pypospack.task.lammps import LammpsSimulationError

#---- imports for IterativeSampler
from pypospack.pyposmat.data.logfile import PyposmatLogFile

# --- import for PyposmatEngine
from pypospack.pyposmat.engines.mc_sampler import PyposmatBadParameterError


class PyposmatEngine(object):
    """
        Args:
            filename_in(str):
            filename_out(str):
            base_directory(str): This is the base directory from which the
                PyposmatEngine will create and run simulations.  By default
                this is set to None, which means it will use the current
                working directory as the base directory.
            fullauto(bool):
        Attributes:
            pyposmat_filename_in(str)
            pyposmat_filename_out(str)
            base_directory(str)
            rank_directory(str): This reflect the MPI rank of the processsor
                that the PyposmatEngine is running on.  If there is no MPI
                available, this is automatically set to rank0000.
            configuration(pypospack.pyposmat.PyposmatConfigurationFile)
            qoi_manager(pypospack.qoi.QoiManager)
            task_mamanger(pypospack.task.TaskManager)
    """
    def __init__(self,
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.results.out',
            base_directory = None,
            fullauto = False):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        self.pyposmat_filename_in = filename_in
        self.pyposmat_filename_out = filename_out

        self.base_directory = None
        self.rank_directory = None
        self.configuration = None
        self.qoi_manager = None
        self.task_manager = None

        if base_directory is None:
            self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
        else:
            msg_err = "base_directory has to be a string"
            raise ValueError(msg_err)

        if fullauto:
            self.configure()

    @property
    def structures(self):
        """(collections.OrderedDict)"""
        return self.configuration.structures

    @property
    def potential(self):
        """(collections.OrderedDict)"""
        return self.configuration.potential

    def configure(self):
        """

        When writing a new PypospackEngine this method will likely have
        to be modified
        """
        self.create_base_directories()
        self.read_configuration_file()
        self.configure_qoi_manager()
        self.configure_task_manager()

    def create_base_directories(self,base_directory=None):
        assert isinstance(base_directory,str) or base_directory is None

        # <-------- determine the base directory.
        if base_directory is None:
            if self.base_directory is None:
                self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
        else:
            msg_err = "the base directory must be a string"
            raise ValueError(msg_err)

        # <-------- create the base directory if the base directory does
        #           not exist
        if not os.path.exists(self.base_directory):
            os.mkdirs(self.base_directory)

        # <-------- the rank directory is determined by the MPI rank
        #           this is not implemented yet
        if self.rank_directory is None:
            _rank_directory = "rank0"
            self.rank_directory = os.path.join(
                    self.base_directory,
                    _rank_directory)


    def read_configuration_file(self,filename=None):
        assert isinstance(filename,str) or filename is None

        _filename_in = None
        if filename is None:
            _filename_in = self.pyposmat_filename_in
        else:
            _filename_in = filename
            self.pyposmat_filename_in = filename

        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(_filename_in)

    def configure_qoi_manager(self,qois=None):
        if qois is None:
            _qois= self.configuration.qois

        self.qoi_manager = QoiManager(qoi_database=_qois,fullauto=True)


    def configure_task_manager(self):
        # <-------- local variables
        _base_directory = self.base_directory
        _tasks = self.qoi_manager.tasks
        _structures = self.structures

        # <-------- configure task manager
        self.task_manager = TaskManager(
                base_directory=_base_directory)
        self.task_manager.configure(
                tasks = _tasks,
                structures = _structures)

    def evaluate_parameter_set(self,parameters):
        self.configure_task_manager()
        _parameters = copy.deepcopy(parameters)
        _potential = copy.deepcopy(self.configuration.potential)
        
        try:
            self.task_manager.evaluate_tasks(
                    parameters=_parameters,
                    potential=_potential)
        except LammpsSimulationError as e:
            str_neighlist_overflow = 'Neighbor list overflow'
            raise
        except:
            print("--- FATAL ERROR ---")
            print("self.configuration.potential:")
            for k,v in self.configuration.potential.items():
                print("\t",k,'=',v)
            print("current_parameter_set:")
            for k,v in _parameters.items():
                print("\t",k,'=',v)
            print("--- END ERROR INFO ---")

            print(type(self.configuration.potential))
            raise
        else:
            # send the results from the task calculations to calculate QOIs
            _task_results = self.task_manager.results
            self.qoi_manager.calculate_qois(
                    task_results=_task_results)

            # populate qoi values
            _qoi_results = OrderedDict()
            for k_qoi,v_qoi in self.qoi_manager.qois.items():
                _qoi_val = v_qoi['qoi_val']
                _qoi_results[k_qoi] = _qoi_val

            # populate errors
            _qoi_errors = OrderedDict()
            for k_qoi,v_qoi in self.qoi_manager.qois.items():
                _qoi_error_name = '{}.{}'.format(k_qoi,'err')
                _qoi_error = v_qoi['qoi_err']
                _qoi_errors[_qoi_error_name] = _qoi_error

            _results = OrderedDict()
            _results['parameters'] = copy.deepcopy(_parameters)
            _results['qois'] = copy.deepcopy(_qoi_results)
            _results['errors'] = copy.deepcopy(_qoi_errors)

        return _results


class PyposmatClusterSampler(PyposmatEngine):

    def __init__(self, log,
            filename_in='pyposmat.config.in',
            filename_out='pyposmat.results.out',
            mpi_rank=None,
            mpi_size=None,
            base_directory=None):

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)
        self.mpi_rank=mpi_rank
        self.mpi_size=mpi_size
        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out
        self.pyposmat_data_bad_filename = 'pypospack.results.bad'
        # for cluster analysis
        self.pyposmat_configuration_fn = filename_in
        if base_directory is None:
            self.base_directory = os.getcwd()
        self.log = log

    @property
    def configuration_fn(self): return self.pyposmat_configuration_fn

    @property
    def data_fn(self): return self.pyposmat_data_in_filename

    def configure_pyposmat_datafile_in(self,filename):
        self.pyposmat_data_in_filename = filename
        self.pyposmat_datafile_in = PyposmatDataFile(filename)
        self.data = PyposmatDataFile()
        self.data.read(filename)

    def configure_pyposmat_datafile_out(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.pyposmat_data_out_filename = filename
        self.pyposmat_datafile_out = PyposmatDataFile(filename)

    def read_configuration_file(self,filename=None):
        if filename is not None:
            _filename = filename
            self.pyposmat_configuration_fn = filename
        else:
            _filename = self.pyposmat_configuration_fn

        try:
            PyposmatEngine.read_configuration_file(self,filename=_filename)
        except FileNotFoundError as e:
            print("Cannot read filename:")
            print("    filename={}".format(filename))
            print("    o.pyposmat_configuration_fn={}".format(self.pyposmat_configuration_fn))
            print("    _filename={}".format(_filename))
            raise e

        # get information from the PyposmatConfigurationFile object
        self.structure_directory = self.configuration.structures['structure_directory']
        self.parameter_names = [p for p in self.configuration.sampling_distribution]
        self.qoi_names = [k for k in self.configuration.qois]
        self.error_names = ['{}.err'.format(k) for k in self.qoi_names]
        self.parameter_distribution_definition =\
                self.configuration.sampling_distribution
        
        try:
            self.free_parameter_names = [k for k,v in self.parameter_distribution_definition.items() if v[0] != 'equals']
        except KeyError as e:
            print(self.parameter_distribution_definition.items())
            raise

        if self.configuration.sampling_constraints is not None:
            self.parameter_constraints = copy.deepcopy(self.configuration.sampling_constraints)
        else:
            self.parameter_constraints = OrderedDict()

        self.constrained_parameter_names = []
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                self.constrained_parameter_names.append(p)

    def run_simulations(self,
            i_iteration,
            n_samples=None,
            filename=None):
        # process the arguments of the method first

        # arg: i_iteration
        i = i_iteration

        self.log.write("Processing n_samples arg in PyposmatClusterSampler...")
        # arg: n_samples
        if n_samples is not None:
            _n_samples = n_samples
        else:
            try:
               _n_samples = self.configuration.sampling_type[i]['n_samples_per_cluster']
            except KeyError as e:
                print("must use \"n_samples_per_cluster\" keyword to describe the number of simulations per cluster")
                raise e

        self.log.write("Processing filename arg in PyposmatClusterSampler...")
        # arg: filename
        if filename is None:
            _filename = self.pyposmat_data_in_filename
        else:
            _filename = filename

        self.log.write('Reading the data file in PyposmatClusterSampler...')
        # read in the the datafile
        self.data = PyposmatDataFile()
        self.data.read(_filename)

        # determine the sampling type
        _sampling_type = self.configuration.sampling_type[i]['type']
        if _sampling_type == 'kde_w_clusters':
            if 'cluster_id' in self.data.df.columns:
                print('datafile already has clusters...')
            else:
                print('clustering data from:\n{}'.format(_filename))
                # prepare the clustering params in an ordered dict
                
                cluster_args = None
                print(self.configuration.sampling_type[i])
                if 'cluster_args' not in self.configuration.sampling_type[i]:
                    cluster_args = self.get_clustering_parameters(
                            configuration_fn=self.configuration_fn,
                            data_fn=_filename)
                else:
                    cluster_args = self.configuration.sampling_type[i]['cluster_args']
                    cluster_args['configuration_fn'] = os.path.join(os.getcwd(), self.configuration_fn)
                    cluster_args['data_fn'] = os.path.join(os.getcwd(), _filename)
                # send the data to be clustered
                obj_cluster_analysis = PyposmatClusterAnalysis.init_from_ordered_dict(cluster_args)
                obj_cluster_analysis.preprocess_data(cluster_args)
                obj_cluster_analysis.calculate_manifold(cluster_args)
                obj_cluster_analysis.calculate_kNN_analysis(cluster_args)
                obj_cluster_analysis.calculate_clusters(cluster_args)
                # BETTER IMPLEMENTATION FOR CHECKING DECOMPOSITION FAILURE
                # obj_cluster_analysis.isValidPartition does not yet work
                '''
                while True:
                    if not obj_cluster_analysis.isValidPartition():
                        obj_cluster_analysis.calculate_manifold(cluster_args)
                        obj_cluster_analysis.calculate_kNN_analysis(cluster_args)
                        obj_cluster_analysis.calculate_clusters(cluster_args)
                    else:
                        break
                '''
                # use newly clustered data for sampling
                self.data.df = obj_cluster_analysis.data.df
        else:
            raise ValueError(
                'unknown sampling type:{}'.format(
                    _sampling_type
                )
            )
        # actual usage of the clustered data goes here
        # get unique cluster ids
        while True:
            try:
                self._run_mc_cluster_sampling(i=i,
                                              _n_samples=_n_samples)
            except np.linalg.linalg.LinAlgError:
                self.log.write("Encountered LinAlgError in _run_mc_cluster_sampler. Rebuilding partition...")
                if 'cluster_args' not in self.configuration.sampling_type[i]:
                    cluster_args = self.get_clustering_parameters(
                            configuration_fn=self.configuration_fn,
                            data_fn=_filename)
                else:
                    cluster_args = self.configuration.sampling_type[i]['cluster_args']
                    cluster_args['configuration_fn'] = os.path.join(os.getcwd(), self.configuration_fn)
                    cluster_args['data_fn'] = os.path.join(os.getcwd(), _filename)
                obj_cluster_analysis = PyposmatClusterAnalysis.init_from_ordered_dict(cluster_args)
                obj_cluster_analysis.preprocess_data(cluster_args)
                obj_cluster_analysis.calculate_manifold(cluster_args)
                obj_cluster_analysis.calculate_kNN_analysis(cluster_args)
                obj_cluster_analysis.calculate_clusters(cluster_args)
                self.data.df = obj_cluster_analysis.data.df
            else:
                break

        if self.mpi_rank == 0:
            print(i_iteration,_n_samples,_sampling_type)

    def _run_mc_cluster_sampling(self, i, _n_samples):
        cluster_ids = set(self.data.df['cluster_id'])
        self.log.write("cluster_ids={}".format(cluster_ids))
        for cluster_id in cluster_ids:
            self.log.write("cluster_id={}".format(cluster_id))
            mc_sampler = PyposmatMonteCarloSampler(log=self.log)
            # mc_sampler.pyposmat_filename_out = 'iammadeofmagic.out'
            mc_sampler.configuration = PyposmatConfigurationFile()
            mc_sampler.configuration.read(self.configuration_fn)
            mc_sampler.configuration.sampling_type[i] = OrderedDict()
            mc_sampler.configuration.sampling_type[i]['type'] = 'kde'
            mc_sampler.configuration.sampling_type[i]['n_samples'] = _n_samples

            mc_sampler.create_base_directories()
            mc_sampler.read_configuration_file(self.configuration_fn)
            _structure_dir = mc_sampler.configuration.structures['structure_directory']
            mc_sampler.configuration.structures['structure_directory'] = \
                os.path.join('..', _structure_dir)
            mc_sampler.configure_qoi_manager()
            mc_sampler.configure_task_manager()
            mc_sampler.configure_pyposmat_datafile_out()
            # pyposmat_datafile_out = PyposmatDataFile(filename_out)

            mc_sampler.print_structure_database()
            mc_sampler.print_sampling_configuration()
            mc_sampler.print_initial_parameter_distribution()

            mc_sampler.run_kde_sampling(n_samples=_n_samples,
                                        filename_in=self.data_fn,
                                        cluster_id=cluster_id)
            self.log.write("Finished mc_sampler.run_kde_sampling of cluster_id={}".format(cluster_id))

    def get_clustering_parameters(
            self,
            configuration_fn,
            data_fn
            ):

        d = OrderedDict()
        d['configuration_fn'] = configuration_fn
        d['data_fn'] = data_fn
        d['include_parameters'] = True
        d['include_qois'] = True
        d['include_errors'] = False

        d['preprocessing'] = OrderedDict()
        d['preprocessing']['type'] = 'standard_scaler'
        d['preprocessing']['args'] = OrderedDict()
        d['preprocessing']['args']['copy'] = True
        d['preprocessing']['args']['with_mean'] = True
        d['preprocessing']['args']['with_std'] = True

        d['manifold'] = OrderedDict()
        d['manifold']['type'] = 'tsne'
        d['manifold']['args'] = OrderedDict()
        d['manifold']['args']['n_components'] = 2
        d['manifold']['args']['perplexity'] = 30
        d['manifold']['args']['early_exaggeration'] = 12
        d['manifold']['args']['learning_rate'] = 200
        d['manifold']['args']['n_iter'] = 5000
        d['manifold']['args']['n_iter_without_progress'] = 300,
        d['manifold']['args']['min_grad_norm'] = 1e-7,
        # d['manifold']['args']['metric']='euclidean',
        d['manifold']['args']['init'] = 'pca',
        d['manifold']['args']['verbose'] = 0,
        d['manifold']['args']['random_state'] = None
        # method='barnes_hut'
        # angle=0.5

        d['neighbors'] = OrderedDict()
        d['neighbors']['type'] = 'ball_tree'
        d['neighbors']['kNN'] = 4
        d['neighbors']['args'] = OrderedDict()
        d['neighbors']['args']['leaf_size'] = 40
        d['neighbors']['args']['metric'] = 'minkowski'

        d['cluster'] = OrderedDict()
        d['cluster']['type'] = 'dbscan'
        d['cluster']['args'] = OrderedDict()
        d['cluster']['args']['eps'] = OrderedDict()
        d['cluster']['args']['eps']['NN'] = 3
        d['cluster']['args']['eps']['percentile'] = .99
        d['cluster']['args']['min_samples'] = 10
        d['cluster']['args']['metric'] = 'euclidean'
        d['cluster']['args']['metric_params'] = None
        d['cluster']['args']['algorithm'] = 'auto'
        d['cluster']['args']['leaf_size'] = 30
        d['cluster']['args']['p'] = None
        # add conditional here to check for custom clustering params in the config
        # pseudo: if 'clustering_params' in self.configuration:
        # pseudo:    for k,v in self.configuration['clustering_params'].items():
        # pseudo:        d[k] = v
        # not yet implemented in configuration
        return d


class PyposmatMonteCarloSampler(PyposmatEngine):
    def __init__(self, log,
            filename_in='pypospack.config.in',
            filename_out='pypospack.results.out',
            mpi_rank=None,
            mpi_size=None,
            base_directory=None):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        assert type(base_directory) in [str,type(None)]

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)
        self.mpi_rank=mpi_rank
        self.mpi_size=mpi_size
        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out
        self.pyposmat_data_bad_filename = 'pypospack.results.bad'
        self.log = log

    def _log(self,str_msg):
        print(str_msg)

    def configure_pyposmat_datafile_in(self,filename):
        self.pyposmat_data_in_filename = filename
        self.pyposmat_datafile_in = PyposmatDataFile(filename)

    def configure_pyposmat_datafile_out(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.pyposmat_data_out_filename = filename
        self.pyposmat_datafile_out = PyposmatDataFile(filename)

    def read_configuration_file(self,filename=None):
        PyposmatEngine.read_configuration_file(self,filename=filename)
        self.structure_directory = self.configuration.structures['structure_directory']
        self.n_iterations = self.configuration.sampling_type['n_iterations']
        self.parameter_names = [p for p in self.configuration.sampling_distribution]
        self.qoi_names = [k for k in self.configuration.qois]
        self.error_names = ['{}.err'.format(k) for k in self.qoi_names]
        self.parameter_distribution_definition =\
                self.configuration.sampling_distribution
        
        try:
            self.free_parameter_names = [k for k,v in self.parameter_distribution_definition.items() if v[0] != 'equals']
        except KeyError as e:
            print(self.parameter_distribution_definition.items())
            raise
        if self.configuration.sampling_constraints is not None:
            self.parameter_constraints = copy.deepcopy(self.configuration.sampling_constraints)
        else:
            self.parameter_constraints = OrderedDict()

        self.constrained_parameter_names = []
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                self.constrained_parameter_names.append(p)

    def run_simulations(self,i_iteration,n_samples=None,filename=None):
        i = i_iteration
        _sampling_type = self.configuration.sampling_type[i]['type']
        _n_samples = self.configuration.sampling_type[i]['n_samples']
        
        if n_samples is not None:
            _n_samples = n_samples

        if _sampling_type == 'parametric':
            self.run_parameteric_sampling(n_samples=_n_samples)
        elif _sampling_type == 'kde':
            if filename is None:
                raise ValueError('cannot do kde sampling with out filename')
            self.run_kde_sampling(n_samples=_n_samples,filename_in=filename)
        elif _sampling_type == 'from_file':
            if filename is None:
                raise ValueError('cannot do filesampling without file')
            self.run_file_sampling(filename)
        else:
            raise ValueError(
                'unknown sampling type:{}'.format(
                    _sampling_type
                )
            )
        if self.mpi_rank == 0:
            print(i_iteration,_n_samples,_sampling_type)

    def run_parameteric_sampling(self,n_samples):
        _rv_generators = OrderedDict()
        for p in self.free_parameter_names:
            distribution_type = self.parameter_distribution_definition[p][0]
            if distribution_type == 'uniform':
                _a = self.parameter_distribution_definition[p][1]['a']
                _b = self.parameter_distribution_definition[p][1]['b']
                _loc = _a
                _scale = _b-_a
                _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
            else:
                raise ValueError('unknown distribution type: {}'.format(
                    distribution_type))

        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration = time.time()
        _n_errors = 0

        for i_sample in range(n_samples):
            # generate parameter set
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])

            # generate free parameters
            for p in self.free_parameter_names:
                _parameters[p] = _rv_generators[p].rvs(size=1)[0]

            # generate constrained parameters
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
                        _parameters[p] = eval(_str_eval)

            # generate wierd things
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        # required for EAM potentials to calculate dens_max for embedding function
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)

            try:
                # check constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))

                    try:
                        _is_constraint_ok = eval(_eval_str)
                    except NameError as e:
                        _str = str(e)
                        _regex_str = "name \'d_NiNi_r0\' is not defined"
                        _err_msg = "BadQoiConstraint:\n"
                        _err_msg += "\t{}".format(k)
                        _err_msg += "\t{}".format(v)
                        _err_msg += "\t{}".format(_eval_str)
                        self._log(_err_msg)
                        raise
                    if eval(_eval_str) is False:
                        raise PyposmatBadParameterError()

                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print(_str_msg)
                    sys.stdout.flush()

    def run_kde_sampling(self,n_samples,filename_in,cluster_id=None):

        _datafile_in = None
        if cluster_id is None:
            _datafile_in = PyposmatDataFile()
            _datafile_in.read(filename_in)
        else:
            _datafile_in = PyposmatDataFile()
            _datafile_in.read(filename_in)

            # redefining the dataframe by subselecting the cluster_id we are interested in
            _datafile_in.df = _datafile_in.df.loc[
                _datafile_in.df['cluster_id'] == cluster_id
            ]

        _X = _datafile_in.df[self.free_parameter_names].loc[_datafile_in.df['cluster_id'] == cluster_id].values.T

        try:
            _h = Chiu1999_h(_X)
            kde_bw_type = 'Chiu1999'
        except LinAlgError as e:
            print('filename:{}'.format(filename_in))
            raise

        d = OrderedDict()
        d['kde_bandwidth'] = OrderedDict()
        d['kde_bandwidth']['type'] = kde_bw_type
        d['kde_bandwidth']['h'] = _h

        _rv_generator = scipy.stats.gaussian_kde(_X,_h)
        print('Chiu1999_h:{}'.format(_h))

        self.pyposmat_datafile_out.df = copy.deepcopy(_datafile_in.df)
        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration = time.time()
        _n_errors = 0
        self.log.write("Begin iteration over n_samples in PyposmatMonteCarloSampler...")
        for i_sample in range(n_samples):
            print(i_sample)
            # generate parameter set
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])
            _free_parameters = _rv_generator.resample(1)
            for i,v in enumerate(self.free_parameter_names):
                _parameters[v] = float(_free_parameters[i,0])

            # generate constrained parameters
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
                        _parameters[p] = eval(_str_eval)

            # generate wierd things
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)
            try:
                # check constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))
                    if eval(_eval_str) is False:
                        raise PyposmatBadParameterError()

                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print(_str_msg)
                    self.log.write(_str_msg)

    def run_file_sampling(self,filename_in):

        _datafile_in = PyposmatDataFile(filename=filename_in)
        _datafile_in.read()
        # configure random number generator

        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration = time.time()
        if self.mpi_rank is None:
            self.mpi_rank = 0
        if self.mpi_size is None:
            self.mpi_size = 1

        _n_errors = 0
        i_sample = 0
        for row in _datafile_in.df.iterrows():
            if self.mpi_rank != row[0]%self.mpi_size:
                continue
            _parameters = OrderedDict([(p,row[1][p]) for p in self.parameter_names])

            # generate wierd things
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)
            try:
                # check constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))
                    if eval(_eval_str) is False:
                        raise PyposmatBadParameterError()

                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                i_sample = i_sample+1
                if (i_sample)%10 == 0:
                    n_samples_completed = i_sample
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print('rank{}:'.format(self.mpi_rank)+_str_msg)

    def calculate_equilibrium_density(self,a0,latt,parameters):
        _parameters = OrderedDict()
        for k,v in parameters.items():
            if k.startswith('d_'):
                _parameters[k[2:]] = v
            s = k[2:].split('_')[0]
        _potential_type = self.configuration.potential['density_type']
        _symbols = self.configuration.potential['symbols']
        _module_name,_class_name = PotentialObjectMap(
                potential_type=_potential_type)
        try:
            _module = importlib.import_module(_module_name)
            _class = getattr(_module,_class_name)
            _dens_potential = _class(symbols=_symbols)
        except:
            raise

        if latt == 'fcc':
            d = OrderedDict([
                ('1NN',2/(2**0.5)*a0),
                ('2NN',1.000*a0),
                ('3NN',1.225*a0)])
            Z= OrderedDict([
                ('1NN',12),
                ('2NN',6),
                ('3NN',24)])
            rcut = (d['2NN']+d['3NN'])/2.

            rmax = 10.
            r = np.linspace(1,10,5000)*rmax/10
            rho = _dens_potential.evaluate(r,_parameters,rcut)

            rho_e = 0
            for m in Z:
                if d[m] < rcut:
                    rho_e += Z[m]*np.interp(d[m],r,rho[s])

            return rho_e

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
            _sample_type = self.configuration.sampling_type[i]['type']
            if _sample_type == 'kde_w_clusters':
                _n_samples = self.configuration.sampling_type[i]['n_samples_per_cluster']
            else:
                _n_samples = self.configuration.sampling_type[i]['n_samples']
            print('{:^10} {:^10} {:^20}'.format(i,_n_samples,_sample_type))

    def print_initial_parameter_distribution(self):
        print(80*'-')
        print('{:80}'.format('INITIAL PARAMETER DISTRIBUTION'))
        print(80*'-')
        for p in self.parameter_distribution_definition:
            if p in self.free_parameter_names:
                str_free = 'free'
                print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
                    p,
                    str_free,
                    self.parameter_distribution_definition[p][0],
                    self.parameter_distribution_definition[p][1]['a'],
                    self.parameter_distribution_definition[p][1]['b']))
            else:
                str_free = 'not_free'
                print('{:^20} {:^10}'.format(p,str_free))


class PyposmatIterativeSampler(object):
    def __init__(self,
            configuration_filename,is_restart=False):
        self.RANK_DIR_FORMAT = 'rank_{}'
        self.mpi_comm = None
        self.mpi_rank = None
        self.mpi_size = None
        self.mpi_nprocs = None
        self.n_iterations = None
        self.rv_seed = None
        self.rv_seeds = None

        self.configuration_filename = configuration_filename
        self.configuration = None
        self.mc_sampler = None

        self.root_directory = os.getcwd()
        self.data_directory = 'data'
        self.is_restart = is_restart
        self.start_iteration = 0

        self.log_fn = os.path.join(self.root_directory, self.data_directory, 'pyposmat.log')
        self.log = PyposmatLogFile(filename=self.log_fn)

    def run_restart(self):
        if self.configuration is None:
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(self.configuration_filename)


        # determine if there was a seed file
        _init_fn = self.find_initial_parameters_files()

        # get contents of the data directory if it exists
        _data_dir = os.path.join(self.root_directory,self.data_directory)
        self.i_iterations, _data_dir_fns = self.analyze_rank_directories(
                data_dir=_data_dir)

        # get contents of the rank directories
        _root_dir = self.root_directory
        n_ranks, _rank_data_fns = self.analyze_rank_directories(
                root_dir=_root_dir)

        if self.mpi_rank == 0:
            pass
        pass

    def run_all(self):
        self.setup_mpi_environment()
        self.determine_rv_seeds()
        MPI.COMM_WORLD.Barrier()

        self.start_iteration = 0
        for i in range(self.start_iteration,self.n_iterations):
            if self.mpi_rank == 0:
                self.log.write(80*'-')
                self.log.write('{:80}'.format('BEGIN ITERATION {}/{}'.format(
                               i+1, self.n_iterations)))
                self.log.write(80*'-')
            MPI.COMM_WORLD.Barrier()

            self.run_simulations(i)
            MPI.COMM_WORLD.Barrier()
            self.log.write('simulations complete...')

            if self.mpi_rank == 0:
                self.log.write('merging files...')
                self.merge_files(i)
                self.log.write('analyzing results...')
                self.analyze_results(i)
            MPI.COMM_WORLD.Barrier()
        self.log.write(80*'-')
        self.log.write('JOBCOMPLETE')

    def run_simulations(self,i_iteration):
        self.rank_directory = self.RANK_DIR_FORMAT.format(
                self.mpi_rank)

        # if the directory exists delete it
        if os.path.isdir(self.rank_directory):
            shutil.rmtree(self.rank_directory)
        os.makedirs(self.rank_directory)

        # change execution context for this rank
        # this provides a directory for each worker directory so that the
        # disk IO writes don't conflict
        os.chdir(self.rank_directory)

        _config_filename = os.path.join(
                self.root_directory,
                self.configuration_filename)

        _results_filename = os.path.join(
                self.root_directory,
                self.rank_directory,
                'pyposmat.results.out')

        # set random seed
        np.random.seed(self.rv_seeds[self.mpi_rank,i_iteration])

        self.log.write("Initializing PyposmatMonteCarloSampler...")
        # initialize()
        self.mc_sampler = PyposmatMonteCarloSampler(
                filename_in = _config_filename,
                filename_out = _results_filename,
                mpi_rank = self.mpi_rank,
                mpi_size = self.mpi_size,
                log=self.log)
        self.mc_sampler.create_base_directories()
        self.mc_sampler.read_configuration_file()
        _structure_dir = self.mc_sampler.configuration.structures['structure_directory']
        self.mc_sampler.configuration.structures['structure_directory'] = \
                os.path.join('..',_structure_dir)
        self.mc_sampler.configure_qoi_manager()
        self.mc_sampler.configure_task_manager()
        self.mc_sampler.configure_pyposmat_datafile_out()
        #pyposmat_datafile_out = PyposmatDataFile(filename_out)

        if self.mpi_rank == 0:
            self.mc_sampler.print_structure_database()
            self.mc_sampler.print_sampling_configuration()
        if self.mpi_rank == 0 and i_iteration == 0:
            self.mc_sampler.print_initial_parameter_distribution()
        if self.mpi_rank == 0:
            self.log.write(80*'-')
        MPI.COMM_WORLD.Barrier()

        _mc_config = self.mc_sampler.configuration.sampling_type[i_iteration]
        # choose sampling type
        _mc_sample_type = _mc_config['type']
        self.log.write("_mc_sample_type={}".format(_mc_sample_type))
        
        # <----- paramter sampling type ---------------------------------------
        if _mc_sample_type == 'parametric':
            _mc_n_samples = _mc_config['n_samples']
            # determine number of sims for this rank
            _n_samples_per_rank = int(_mc_n_samples/self.mpi_size)
            if _mc_n_samples%self.mpi_size > self.mpi_rank:
                _n_samples_per_rank += 1
            self.mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank)
        
        # <----- kde sampling sampling type ---------------------------------------
        elif _mc_sample_type == 'kde':
            _mc_n_samples = _mc_config['n_samples']
            # determine number of sims for this rank
            _n_samples_per_rank = int(_mc_n_samples/self.mpi_size)
            if _mc_n_samples%self.mpi_size > self.mpi_rank:
                _n_samples_per_rank += 1
            _filename_in = ''
            if 'file' in _mc_config:
                _filename_in = os.path.join(
                    self.root_directory,
                    _mc_config['file']
                )
            else:
                _filename_in = os.path.join(
                    self.root_directory,
                    self.data_directory,
                    'pyposmat.kde.{}.out'.format(i_iteration))

            self.mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank,
                    filename=_filename_in)

        # <----- sampling from a file type ---------------------------------------
        # get parameters from file
        elif _mc_sample_type == 'from_file':
            _mc_n_samples = _mc_config['n_samples']
            
            # determine number of sims for this rank
            _n_samples_per_rank = int(_mc_n_samples/self.mpi_size)
            if _mc_n_samples%self.mpi_size > self.mpi_rank:
                _n_samples_per_rank += 1
            
            _filename_in = os.path.join(
                self.root_directory,
                _mc_config['file']
            )
            
            self.mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank,
                    filename=_filename_in
            )
        # <----- kde with clusters sampling type ---------------------------------------
        elif _mc_sample_type == 'kde_w_clusters':
            
            pyposmat_datafile_in = os.path.join(
                self.root_directory,
                self.data_directory,
                "pyposmat.cluster.0.out"
            )

            _config_filename = os.path.join(
                self.root_directory,
                self.configuration_filename)

            # determine number of sims for this rank
            _mc_n_samples = _mc_config['n_samples_per_cluster']
            _n_samples_per_rank = int(_mc_n_samples / self.mpi_size)
            if _mc_n_samples % self.mpi_size > self.mpi_rank:
                _n_samples_per_rank += 1

            self.log.write("Initializing PyposmatClusterSampler...")
            # initialize sampling object
            o = PyposmatClusterSampler(log=self.log)
            o.create_base_directories()
            o.read_configuration_file(filename=_config_filename)
            o.configure_pyposmat_datafile_in(filename=pyposmat_datafile_in)
            
            # fix relative path to structure databae folder
            _structure_dir = o.configuration.structures['structure_directory']
            o.configuration.structures['structure_directory'] = \
                    os.path.join('..',_structure_dir)

            # finish the rest of the initialization
            o.configure_qoi_manager()
            o.configure_task_manager()
            o.configure_pyposmat_datafile_out()

            # run simulations
            self.log.write("Running simulations through PyposmatClusterSampler...")
            o.run_simulations(i_iteration=i_iteration,
                              n_samples=_mc_n_samples,
                              filename=pyposmat_datafile_in)
        else:
            m = "unknown sampling type: {}".format(
               _mc_sample_type
            )
            raise ValueError(m)
        
        # return to root directory
        os.chdir(self.root_directory)

    def get_results_dict(self):
        rd = OrderedDict()
        rd['mpi'] = OrderedDict()
        rd['mpi']['size'] = self.mpi_size

    def setup_mpi_environment(self):
        self.mpi_comm = MPI.COMM_WORLD
        self.mpi_rank = self.mpi_comm.Get_rank()
        self.mpi_size = self.mpi_comm.Get_size()
        self.mpi_procname = MPI.Get_processor_name()
        if self.mpi_rank == 0:
            self.print_mpi_environment()

    def print_mpi_environment(self):
        self.log.write(80*'-')
        self.log.write('{:^80}'.format('MPI COMMUNICATION INFORMATION'))
        self.log.write(80 * '-')
        self.log.write('mpi_size={}'.format(self.mpi_size))

    def determine_rv_seeds(self):
        _randint_low = 0
        _randint_high = 2147483647

        # set original seed
        if self.rv_seed is None:
            self.rv_seed = np.random.randint(
                    low=_randint_low,
                    high=_randint_high)
        np.random.seed(self.rv_seed)

        # determine rank seed
        self.rv_seeds = np.random.randint(
            low=0,
            high=2147483647,
            size=(int(self.mpi_size),self.n_iterations)
            )

        if self.mpi_rank == 0:
            self.print_random_seeds()

    def analyze_data_directories(self,data_dir=None):
        _d = data_dir
        i = 0
        contents = []
        if not os.path.exists(_d): return i, contents
        if not os.path.isdir(_d): return i, contents

        while True:
            kde_fn = os.path.join(_d,"pyposmat.kde.{}.out".format(i))
            if os.path.exists(kde_fn):
                contents.append(kde_fn)
            else:
                if i > 0:
                    contents.append(results_fn)
                    break

            results_fn = os.path.join(_d,"pyposmat.results.{}.out".format(i))
            if os.path.exists(results_fn): pass
            else:break
            i = i + 1

        return i, contents

    def analyze_rank_directories(self,root_dir=None):
        i = 0
        contents = []

        if root_dir is None:
            _d = self.root_directory
        else:
            _d = root_directory

        while True:
            rank_dir = os.path.join(_d,"rank_{}".format(i))
            if not os.path.exists(rank_dir): break
            if not os.path.isdir(rank_dir): break
            rank_fn = os.path.join("rank_{}".format(i),"pyposmat.results.out")
            if not os.path.exists(os.path.join(_d,rank_fn)): break
            if not os.path.isfile(os.path.join(_d,rank_fn)):
                break
            else:
                contents.append(rank_fn)
            i = i + 1
        return i, contents

    def find_initial_parameters_file(self):
        if 'file' in self.configuration.sampling_type[0]:
            _init_fn =os.path.join(
                self.root_directory,
                self.configuration.sampling_type[0]['file']
            )
            if os.path.exists(_init_fn):
                if os.path.isfile(_init_fn):
                    return _init_fn
                else:
                    return None

    def merge_pypospack_datafiles(datafile_fns):
        d0 = PyposmatDataFile()
        d0.read(filename=datafile_fns[0])
        df0 = d0.df
        for i in range(1,len(datafile_fns)):
            print("merging {}...".format(datafile_fns[i]))
            d = PyposmatDataFile()
            d.read(filename=datafile_fns[i])
            df = d.df

            df0 = pd.concat([df0,df]).drop_duplicates().reset_index(drop=True)
        d0.df = df0
        return d0

    def merge_files(self,i_iteration):

        _dir = self.data_directory
        _n_ranks = self.mpi_size

        datafile = None
        # filename of old kde file
        _filename_kde = os.path.join(
            _dir,'pyposmat.kde.{}.out'.format(i_iteration))
        
        self.log.write('Looking for previous kde file')
        self.log.write('    {}'.format(_filename_kde))
        
        datafile_fns = []
        if os.path.exists(_filename_kde):
            if os.path.isfile(_filename_kde):
                datafile_fns.append(_filename_kde)
        for i_rank in range(_n_ranks):
            rank_fn = os.path.join(
                'rank_{}'.format(i_rank),
                'pyposmat.results.out')
            datafile_fns.append(rank_fn)
        
        names = ['sim_id']\
                + self.parameter_names\
                + self.qoi_names\
                + self.error_names
        types = ['sim_id']\
                + ['param']*len(self.parameter_names)\
                + ['qoi']*len(self.qoi_names)\
                + ['err']*len(self.error_names)
        
        dataframes = OrderedDict()
        for fn in datafile_fns:
            datafile = PyposmatDataFile()
            datafile.read(fn)
            #if fn.startswith('rank') 
            #datafile.df['sim_id'] = datafile.df.apply(
            #    lambda x:"{}_{}_{}".format(
            #        i_iteration,i_rank,str(x['sim_id'])))
            dataframes[fn] = datafile.df[names]
        
              
        df = pd.concat(dataframes).reset_index(drop=True)
        datafile = PyposmatDataFile()
        datafile.df = df
        datafile.parameter_names = self.parameter_names
        datafile.error_names = self.error_names
        datafile.qoi_names = self.qoi_names
        datafile.names = names
        datafile.types = types
        try:
            fn_out = os.path.join(
                _dir,'pyposmat.results.{}.out'.format(i_iteration))
            datafile.write(filename=fn_out)
        except FileNotFoundError as e:
            if not os.path.exists(self.data_directory):
                os.mkdir(self.data_directory)
                datafile.write(filename_fn_out)
            else: raise


    def analyze_results(self,i_iteration):
        data_fn = os.path.join(\
                self.root_directory,
                self.data_directory,
                'pyposmat.results.{}.out'.format(i_iteration))
        config_fn = os.path.join(\
                self.root_directory,
                self.configuration_filename)
        kde_fn = os.path.join(\
                self.root_directory,
                self.data_directory,
                'pyposmat.kde.{}.out'.format(i_iteration+1))

        data_analyzer = PyposmatDataAnalyzer()
        data_analyzer.read_configuration_file(filename=config_fn)
        data_analyzer.read_data_file(filename=data_fn)
        data_analyzer.write_kde_file(filename=kde_fn)

    def read_configuration_file(self,filename=None):
        assert isinstance(filename,str) or filename is None

        if filename is None:
            _filename_in = self.configuration_filename
        else:
            self.configuration_filename = filename
            _filename_in = filename

        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=_filename_in)
        
        self.n_iterations = self.configuration.n_iterations
        self.qoi_names = self.configuration.qoi_names
        self.error_names = self.configuration.error_names
        self.parameter_names = self.configuration.parameter_names
    
        self.log.write(self.parameter_names)
        self.log.write(self.qoi_names)
        self.log.write(self.error_names)

    def print_random_seeds(self):
            self.log.write(80*'-')
            self.log.write('{:^80}'.format('GENERATED RANDOM SEEDS'))
            self.log.write(80*'-')
            self.log.write('')
            self.log.write('rv_seed={}'.format(self.rv_seed))
            self.log.write('')
            self.log.write('{:^8} {:^8} {:^10}'.format('rank','iter','seed'))
            self.log.write('{} {} {}'.format(8*'-',8*'-',10*'-'))
            for i_rank in range(self.mpi_size):
                for i_iter in range(self.n_iterations):
                    self.log.write('{:^8} {:^8} {:>10}'.format(
                                   i_rank,
                                   i_iter,
                                   self.rv_seeds[i_rank, i_iter]))
