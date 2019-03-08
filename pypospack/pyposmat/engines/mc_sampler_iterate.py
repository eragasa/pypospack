"""

VERSION HISTORY:
    2017 - EJR, original release
    08/2018 = RSU, EJR, added cluster sampling functionality and logging

TODO:
    08/2018 - the pyposmat.results.*.out needs to have the sim_ids cleaned up so we can track information regarding rank_id, cluster_id, and iteration_id.  This is at approximately line 442.  EJR&RSU
"""

import os,shutil,sys
import numpy as np
from mpi4py import MPI
import pandas as pd
from collections import OrderedDict

# --- pyposmat data format imports
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatLogFile

# --- pyposmat analysis algorithms
from pypospack.pyposmat.data import PyposmatDataAnalyzer

# --- pyposmat data samplers
from pypospack.pyposmat.engines import PyposmatMonteCarloSampler
from pypospack.pyposmat.engines import PyposmatFileSampler
from pypospack.pyposmat.engines import PyposmatClusterSampler

class PyposmatIterativeSampler(object):
    """  Iterative Sampler which wraps multiple simulation algorithms.

    This class wraps multiple simulation algorithms so that they can be run in an iterative manner.
    Since this class has so many configuration options, the attributes of this class is set
    by a YAML based configuration file.  The class PyposmatConfigurationFile aids in the creation
    and reading of these options.  These attributes are public and be set programmatically within
    a script.

    Notes:
        config_fn = 'data/pyposmat.config.in'

        engine = PyposmatIterativeSampler(configuration_filename=config_fn)
        engine.read_configuration_file()
        engine.run_all()

    Args:
        configuration_filename(str): the filename of the YAML configuration file
        is_restart(bool,optional): When set to True, this argument controls the restart behavior
            of this class.  By default, is set to False
        is_auto(bool,optional): When set to True, this agument will automatically configure the
            class.  By default this is set to False, mostly because this software is currently in
            development, and this necessary to to write integration testing
        log_fn(str,optional): This the filename path where to set logging, by default it is set
           as `pyposmat.log` contained in the configurable data directory
        log_to_stdout(bool,optional): When set to True, all log messages will be directed to
           standard out as well as the log file

    Attributes:
        mpi_comm(MPI.Intracomm)
        mpi_rank(int)
        mpi_size(int)
        mpi_nprocs(int)
        i_iteration(int)
        n_iterations(int)
        rv_seed(int)
        rv_seeds(np.ndarray)
        configuration_filename = configuration(filename)
        configuration(PyposmatConfigurationFile)
        mc_sampler(PyposmatMonteCarloSampler)
        root_directory(str)
        data_directory(str)
        is_restart(bool)
        start_iteration=0
        
    """

    def __init__(self,
            configuration_filename,
            is_restart=False,
            is_auto=False,
            log_fn=None,
            log_to_stdout=True):

        # formats should not contain a trailing end line chracter
        self.SECTION_HEADER_FORMAT = "\n".join([80*'=',"{:^80}",80*"="])
        self.RANK_DIR_FORMAT = 'rank_{}'


        self.mpi_comm = None
        self.mpi_rank = None
        self.mpi_size = None
        self.mpi_nprocs = None
        self.i_iteration = None
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

        self.log_fn = log_fn
        self.log_to_stdout = log_to_stdout
        self.o_log = None
        self.initialize_logger(log_fn=log_fn,log_to_stdout=log_to_stdout)


    def run_all(self):
        """runs all iterations

        This method runs all iterations

        """
        self.setup_mpi_environment()
        self.initialize_data_directory()
        self.start_iteration = 0
        for i in range(self.start_iteration,self.n_iterations):
            self.i_iteration = i

            # log iteration information
            self.log_iteration_information(i_iteration=i)

            if self.is_restart:
                results_fn = os.path.join(self.data_directory,'pyposmat.results.{}.out'.format(i))
                kde_fn = os.path.join(self.data_directory,'pyposmat.kde.{}.out'.format(i))

                if all([os.path.isfile(results_fn),os.path.isfile(kde_fn)]):
                    self.log('this iterations as already been completed')
                else:
                    break

            self.run_simulations(i)
            MPI.COMM_WORLD.Barrier()

            if self.mpi_rank == 0:
                self.log("ALL SIMULATIONS COMPLETE FOR ALL RANKS")

            if self.mpi_rank == 0:
                self.merge_files(i)
                self.analyze_results(i)
            MPI.COMM_WORLD.Barrier()
        self.log(80*'-')
        self.log('JOBCOMPLETE')

    def initialize_sampler(self,config_fn,results_fn,mpi_rank=None,mpi_size=None,o_log=None):
        """ initialize the sampling object 

        This method initializes the `mc_sampler` attribute with a sampler.

        Note:
            This breakout is part of a larger effort within PYPOSPACK, to have 
            more object-oriented approach for parametric sampling.  The goal 
            eventually is to implement an instance of PyposmatBaseSampler, and 
            allow users of this software library to be able to extend this 
            software by simply extending the base class.
        Args:
            config_fn(str): path to the configuration file
            results_fn(str): path to the results file
            mpi_rank(int,optional): the MPI rank of executing this method
            mpi_size(int,optional): the size of the MPI execution group
            o_log(PyposmatLogFile,str,optional): the log file.  If a string is 
                passed, then the sampling class will initialize a separate log 
                file with the string of path created.  If a log file object is 
                passed, then sampling object will use that instance of the 
                object to log information.  By defaut, it will pass the 
                attribute, `o_log`.
        """
        
        assert type(config_fn) is str
        assert type(results_fn) is str
        assert type(mpi_rank) in [type(None),int]
        assert type(mpi_size) in [type(None),int]
        assert type(o_log) in [type(None),PyposmatLogFile,str]

        # check to see if the paths provided are absolute paths
        assert os.path.isabs(config_fn)
        assert os.path.isabs(results_fn)

        if mpi_rank is None: mpi_rank = self.mpi_rank
        if mpi_size is None: mpi_size = self.mpi_size

        self.mc_sampler = PyposmatMonteCarloSampler(
                filename_in = config_fn,
                filename_out = results_fn,
                mpi_rank = mpi_rank,
                mpi_size = mpi_size,
                o_log=o_log)
        self.mc_sampler.create_base_directories()
        self.mc_sampler.read_configuration_file()
        
        # we have to be able to find the structure directory
        self.mc_sampler.configuration.structures['structure_directory'] = self.structure_directory
        self.mc_sampler.configure_qoi_manager()
        self.mc_sampler.configure_task_manager()
        self.mc_sampler.configure_pyposmat_datafile_out()

        self.log_more_iteration_information()

    def initialize_file_sampler(self,
                                config_fn,
                                results_fn,
                                i_iteration=0,
                                mpi_rank=None,
                                mpi_size=None,
                                o_log=None):
        """ initialize the sampling object 

        This method initializes the `mc_sampler` attribute with a sampler.

        Note:
            This breakout is part of a larger effort within PYPOSPACK, to have 
            more object-oriented approach for parametric sampling.  The goal 
            eventually is to implement an instance of PyposmatBaseSampler, and 
            allow users of this software library to be able to extend this 
            software by simply extending the base class.
        Args:
            config_fn(str): path to the configuration file
            results_fn(str): path to the results file
            i_iteration(int,optional): the iteration to sample the file from,
                by default this is set to zero.
            mpi_rank(int,optional): the MPI rank of executing this method
            mpi_size(int,optional): the size of the MPI execution group
            o_log(PyposmatLogFile,str,optional): the log file.  If a string is 
                passed, then the sampling class will initialize a separate log 
                file with the string of path created.  If a log file object is 
                passed, then sampling object will use that instance of the 
                object to log information.  By defaut, it will pass the 
                attribute, `o_log`.
        """

        assert type(config_fn) is str
        assert type(results_fn) is str
        assert type(mpi_rank) in [type(None),int]
        assert type(mpi_size) in [type(None),int]
        assert type(o_log) in [type(None),PyposmatLogFile,str]

        # check to see if the paths provided are absolute paths
        assert os.path.isabs(config_fn)
        assert os.path.isabs(results_fn)

        if mpi_rank is None: mpi_rank = self.mpi_rank
        if mpi_size is None: mpi_size = self.mpi_size

        # get the absolute path of the datafile we are sampling from
        data_in_fn = None
        if os.path.isabs(self.configuration.sampling_type[i_iteration]['file']):
            data_in_fn = self.configuration.sampling_type[i_iteration]['file']
        else:
            data_in_fn = os.path.join(
                   self.root_directory,
                   self.configuration.sampling_type[i_iteration]['file']
            )

        
        data_out_fn = results_fn

        self.mc_sampler = PyposmatFileSampler(
                config_fn = config_fn,
                data_in_fn = data_in_fn,
                data_out_fn = data_out_fn,
                mpi_rank = mpi_rank,
                mpi_size = mpi_size,
                o_log=o_log,
                fullauto=False)
        self.mc_sampler.create_base_directories()
        self.mc_sampler.read_configuration_file()
        
        # we have to be able to find the structure directory
        self.mc_sampler.configuration.structures['structure_directory'] = self.structure_directory
        self.mc_sampler.configure_qoi_manager()
        self.mc_sampler.configure_task_manager()
        self.mc_sampler.configure_datafile_out()

        self.log_more_iteration_information()

    def initialize_rank_directory(self):
        """ create the rank directory

        This method defines the rank directory as an absolute path and stores it in
        the attribute `rank_directory`.  If a current directory exists there, then
        it is deleted with alll it's contents and then recreated.

        """
        rank_directory = os.path.join(self.root_directory,self.RANK_DIR_FORMAT.format(self.mpi_rank))
      
        # find the directory, delete it and it's constants and then recreates ot
        if os.path.isdir(rank_directory):
            shutil.rmtree(rank_directory)
        os.mkdir(rank_directory)

        self.rank_directory = rank_directory

    def run_simulations(self,i_iteration):
        """ run simulation for a single iteration

        Each rank is given a different execution context so that the disk IO 
        don't conflict
        """
        self.initialize_rank_directory()
        config_filename = self.configuration_filename
        results_filename = os.path.join(self.rank_directory,'pyposmat.results.out')
        bad_parameters_filename = os.path.join(self.rank_directory,'pyposmat.bad_parameters.out')

        # change execution context for this rank
        os.chdir(self.rank_directory)

        # set random seed
        self.determine_rv_seeds()
        self.log_random_seeds(i_iteration=i_iteration)

        sampling_type = self.configuration.sampling_type[i_iteration]['type']
        if self.mpi_rank == 0:
            self.log("sampling_type={}".format(sampling_type))
        
        # <----- parameter sampling type ---------------------------------------
        if sampling_type== 'parametric':
            self.initialize_sampler(
                    config_fn=config_filename,
                    results_fn=results_filename,
                    mpi_rank=self.mpi_rank,
                    mpi_size=self.mpi_size,
                    o_log=self.o_log)

            self.run_parametric_sampling(i_iteration=i_iteration)
        
        # <----- kde sampling sampling type ---------------------------------------
        elif sampling_type == 'kde':
            self.initialize_sampler(
                    config_fn=config_filename,
                    results_fn=results_filename,
                    mpi_rank=self.mpi_rank,
                    mpi_size=self.mpi_size,
                    o_log=self.o_log)
            
            self.run_kde_sampling(i_iteration=i_iteration)

        # <----- sampling from a file type ---------------------------------------
        # get parameters from file
        elif sampling_type == 'from_file':

            self.initialize_file_sampler(
                    config_fn=config_filename,
                    results_fn=results_filename,
                    mpi_rank=self.mpi_rank,
                    mpi_size=self.mpi_size,
                    o_log=self.o_log)
            
            self.run_file_sampling(i_iteration=i_iteration)

        # <----- kde with clusters sampling type ---------------------------------------
        elif _mc_sample_type == 'kde_w_clusters':
            cluster_fn = "pyposmat.cluster.{}.out".format(i_iteration)
            pyposmat_datafile_in = os.path.join(
                self.root_directory,
                self.data_directory,
                cluster_fn
            )
                
            _config_filename = os.path.join(
                self.root_directory,
                self.configuration_filename)

            # determine number of sims for this rank
            _mc_n_samples = _mc_config['n_samples_per_cluster']
            _n_samples_per_rank = int(_mc_n_samples / self.mpi_size)
            if _mc_n_samples % self.mpi_size > self.mpi_rank:
                _n_samples_per_rank += 1

            # initialize sampling object
            o = PyposmatClusterSampler(o_logger=self.log,
                                       mpi_rank=self.mpi_rank,
				       mpi_comm=self.mpi_comm,
				       mpi_size=self.mpi_size)
            o.create_base_directories()
            o.read_configuration_file(filename=_config_filename)
            # check to see if clustered data file exists
            if self.mpi_rank == 0:
                if not os.path.isfile(pyposmat_datafile_in):
                    kde_fn = "pyposmat.kde.{}.out".format(i_iteration)
                    kde_fn = os.path.join(
                            self.root_directory,
                            self.data_directory,
                            kde_fn
                    )
                    o.write_cluster_file(filename=kde_fn, i_iteration=i_iteration)

            MPI.COMM_WORLD.Barrier()

            o.configure_pyposmat_datafile_in(filename=pyposmat_datafile_in)
            # fix relative path to structure databae folder
            _structure_dir = o.configuration.structures['structure_directory']
            o.configuration.structures['structure_directory'] = \
                    os.path.join('..',_structure_dir)
            # finish the rest of the initialization
            o.configure_qoi_manager()
            o.configure_task_manager()
            o.configure_pyposmat_datafile_out()
            MPI.COMM_WORLD.Barrier()
            # run simulations
            o.run_simulations(i_iteration=i_iteration,
                              n_samples=_mc_n_samples,
                              filename=pyposmat_datafile_in)
            MPI.COMM_WORLD.Barrier()
        else:
            m = "unknown sampling type: {}".format(
               _mc_sample_type
            )
            raise ValueError(m)
        
        # return to root directory
        os.chdir(self.root_directory)

    def initialize_data_directory(self,data_directory=None):
        """ determine the absolute path of the data directory and create it

        This method sets the `data_directory` attribute of the class and creates
        the `data directory` if the data directory already exists.

        Args:
            data_directory(str):the path of the data directory, the path can be 
                expressed in either a relative path, or an absolute path
        Returns:
            (str) the absolute path of the data directory
        Raises:
            OSError: if the directory is not able to be created
            
        """

        assert type(data_directory) in [type(None),str]
        assert type(self.data_directory) in [type(None),str]

        # determine the data directory path
        if data_directory is None:
            if self.data_directory is None:
                self.data_directory = os.path.join(self.root_directory,'data')
            else:
                if os.path.isabs(self.data_directory):
                    self.data_directory = data_directory
                else:
                    self.data_directory = os.path.join(self.root_directory,self.data_directory)
        elif os.path.isabs(data_directory):
            # absolute path
            self.data_directory = data_directory
        else:
            # create a absolute path from the relative path
            self.data_directory = os.path.join(self.root_directory,data_directory)
            self.data_directory = os.path.abspath(self.data_directory)

        # create data directory
        try:
            os.mkdir(self.data_directory)
            self.log('created the data directory.')
            self.log('\tdata_directory;{}'.format(self.data_directory))
            return (True, self.data_directory)
        except FileExistsError as e:
            self.log('attempted to create data directory, directory already exists.')
            self.log('\tdata_directory:{}'.format(self.data_directory))
            return (True, self.data_directory)
        except OSError as e:
            self.log('attempted to create data directory, cannot create directory.')
            self.log('\tdata_directory:{}'.format(self.data_directory))
            return (False, self.data_directory)

    def run_parametric_sampling(self,i_iteration):
        """ run parametric sampling 

        Args:
            i_iteration(int): what iteration of the sampling is happening
        """

        assert type(i_iteration) is int
        assert type(self.mc_sampler) is PyposmatMonteCarloSampler

        self.mc_sampler.run_simulations(
                i_iteration=i_iteration,
                n_samples=self.determine_number_of_samples_per_rank(i_iteration=i_iteration))

    def run_kde_sampling(self,i_iteration):
        """ run kde sampling

        Args:
            i_iteration(int): what iteration of the sampling is happening
        """

        assert type(i_iteration) is int
        assert type(self.mc_sampler) is PyposmatMonteCarloSampler

        kde_filename = os.path.join(self.data_directory,
                                    'pyposmat.kde.{}.out'.format(i_iteration))
        n_samples_per_rank = self.determine_number_of_samples_per_rank(i_iteration=i_iteration)

        self.mc_sampler.run_simulations(
                i_iteration=i_iteration,
                n_samples=n_samples_per_rank,
                filename=kde_filename
        )

    def run_file_sampling(self,i_iteration):
        """ run kde sampling

        Args:
            i_iteration(int): the iteration which to sampling for
        """
        assert type(i_iteration) is int
        assert type(self.mc_sampler) is PyposmatFileSampler


        if 'file' in self.configuration.sampling_type[i_iteration]:
            filename = os.path.join(self.root_directory,
                                       self.configuration.sampling_type[i_iteration]['file'])
        else:
            if os.path.isabs(self.data_directory):
                filename = os.path.join(self.data_directory,
                                        'pyposmat.kde.{}.out'.format(i_iteration))
            else:
                filename = os.path,join(self.root_directory,
                                        self.data_directory,
                                        'pyposmat.kde.{}.out'.format(i_iteration))

        if self.mpi_rank == 0:
            self.log(80*'-')
            self.log('{:^80}'.format('file sampling'))
            self.log(80*'-')
            self.log('filename_in:{}'.format(filename))
        MPI.COMM_WORLD.Barrier()

        self.mc_sampler.run_simulations(
                i_iteration=i_iteration,
                n_samples=self.determine_number_of_samples_per_rank(i_iteration=i_iteration),
                filename=filename)


    def determine_number_of_samples_per_rank(self,i_iteration,N_samples=None):
        """ determine the number of samples per rank

        The total number of samples needs to be broken up between the ranks, but roughly
        divided the work evenly.

        Args:
            i_iteration(int): which iteration we are in the simulation
            N_samples(int,optional): the total number of samples we are using for 
                this iteration.  If a number is provided, it will override 
                the number of simulations specified in the configuration file.
        Returns:
            (int): the number of samples for this rank
        """

        assert type(i_iteration) is int
        assert type(N_samples) in [type(None),int]
        assert type(self.configuration) is PyposmatConfigurationFile

        if N_samples is None:
            N_samples = self.configuration.sampling_type[i_iteration]['n_samples']
         
        N_samples_per_rank = int(N_samples/self.mpi_size)
        if N_samples%self.mpi_size > self.mpi_rank:
            N_samples_per_rank += 1

        return N_samples_per_rank

    def initialize_logger(self,log_fn=None,log_to_stdout=None):
        """initialize log object
        
        Args:
            log_fn(str,optional)

        """

        assert type(log_fn) in [type(None),str]
        assert type(log_to_stdout) in [type(None),bool]

        if log_fn is None:
            self.log_fn = os.path.join(self.root_directory, self.data_directory, 'pyposmat.log')
        else:
            self.log_fn = log_fn
        self.o_log = PyposmatLogFile(filename=self.log_fn)
        
        self.log_to_stdout = log_to_stdout

    def setup_mpi_environment(self):
        self.mpi_comm = MPI.COMM_WORLD
        self.mpi_rank = self.mpi_comm.Get_rank()
        self.mpi_size = self.mpi_comm.Get_size()
        self.mpi_procname = MPI.Get_processor_name()
        self.log_mpi_environment()
    
    # random seed management
    def determine_rv_seeds(self,seed=None,i_iteration=None):
        """ set the random variable seed across simulations 
        
        Args:
           seed(int,optional)=a seed to determine the rest of the seeds for
               different ranks and iterations.
        """
        RAND_INT_LOW = 0
        RAND_INT_HIGH = 2147483647

        assert type(seed) in [type(None),int]
        assert type(i_iteration) in [type(None),int]

        if type(i_iteration) is type(None):
            i_iteration = self.i_iteration

        # set the seed attribute
        if type(seed) is int:
            self.rv_seed == seed

        
        # set the seed attribute, if the seed attribute is none
        if self.rv_seed is None:
            self.rv_seed = np.random.randint(
                    low=RAND_INT_LOW,
                    high=RAND_INT_HIGH)
                  
        # if the rv_seed was determined in the script, then all ranks will
        # have the same rv_seed attribute
        np.random.seed(self.rv_seed)

        # each rank, will need it's own seed.  So we sample from the freshly
        # generated random number generator, which is identical across ranks
        self.rv_seeds = np.random.randint(
            low=0,
            high=2147483647,
            size=(int(self.mpi_size),self.n_iterations)
            )

        # now restart the seed for this rank
        np.random.seed(self.rv_seeds[self.mpi_rank,i_iteration])

    # logging methods 
    def log(self,s):
        if self.log_to_stdout: 
            print(s)
        if self.o_log is not None: 
            self.o_log.write(s)

    def log_iteration_information(self,i_iteration):
        """log iteration information
        
        Args:
            i_iteration_id(int):the iteration number
        Returns:
            (str) the log string
        """
        if self.mpi_rank == 0:
            s = self.SECTION_HEADER_FORMAT.format(
                        'Begin Iteration {}/{}'.format(i_iteration+1,self.n_iterations))
            self.log(s)
        MPI.COMM_WORLD.Barrier()
        
        if self.mpi_rank == 0:
            return "\n".join(s)
    
    def log_more_iteration_information(self): 
        #TODO: this logging needs to go into a separate logging method. -EJR
        if self.mpi_rank == 0:
            self.mc_sampler.print_structure_database()
            self.mc_sampler.print_sampling_configuration()
        if self.mpi_rank == 0 and self.i_iteration == 0:
            self.mc_sampler.print_initial_parameter_distribution()
        if self.mpi_rank == 0:
            self.log(80*'-')
        MPI.COMM_WORLD.Barrier()

    def log_mpi_environment(self):
        if self.mpi_rank == 0:
            m = [self.SECTION_HEADER_FORMAT.format('MPI communication information')]
                
            m += ['mpi_size={}'.format(self.mpi_size)]

        MPI.COMM_WORLD.Barrier()

    def log_random_seeds(self,i_iteration):
        if self.mpi_rank == 0:
            self.log(80*'-')
            self.log('{:^80}'.format('GENERATED RANDOM SEEDS'))
            self.log(80*'-')
            self.log('global_seed:{}'.format(str(self.rv_seed)))
            self.log('seeds_for_this_iteration:')
            self.log('{:^8} {:^8}'.format('rank','seed'))
            self.log('{} {}'.format(8*'-',8*'-'))
        MPI.COMM_WORLD.Barrier()
        for i_rank in range(self.mpi_size):
            if self.mpi_rank == i_rank:
                self.log('{:^8} {:>10}'.format(i_rank,self.rv_seeds[i_rank,i_iteration]))
        MPI.COMM_WORLD.Barrier()

    def get_results_dict(self):
        rd = OrderedDict()
        rd['mpi'] = OrderedDict()
        rd['mpi']['size'] = self.mpi_size


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
            if not os.path.exists(rank_dir): 
                break
            if not os.path.isdir(rank_dir): 
                break

            rank_fn = os.path.join("rank_{}".format(i),"pyposmat.results.out")
            if not os.path.exists(os.path.join(_d,rank_fn)):
                break
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
 
    def merge_files(self,i_iteration,last_datafile_fn=None,new_datafile_fn=None):
        """ merge the pyposmat data files

        Args:
            i_iteration(int): the current iteration which just finished
            last_datafile_fn(str,optional): the filename of the last dataset in the data directory.
            new_datafile_fn(str,optional): where to output the file results 
        """

        _dir = self.data_directory

        datafile = None
        
        # filename of old datafile file
        # In default use, 'pyposmat.kde.{i}.out' is used to define the KDE distribution for the ith iteration,
        # Unless the last_datafile_fn is specifically specified, we assume that this is the file used
        if last_datafile_fn is None:
            _datafile_fn_old = os.path.join(_dir,'pyposmat.kde.{}.out'.format(i_iteration))
        else:
            _datafile_fn_old = last_datafile_fn

        # filenames from each rank of the simulations
        # each rank has it's own disk space for the purposes of concurrency, at the end of the iteration it is
        # necessary to concatenate these files
        _n_ranks = self.mpi_size
        _datafiles_fn_ranks = [os.path.join('rank_{}'.format(i),'pyposmat.results.out') for i in range(_n_ranks)]

        # we are concatenating the previous parameterizations and results with the new parameterizations and results
        _datafiles_to_concatenate = []
        if os.path.isfile(_datafile_fn_old):
            _datafiles_to_concatenate.append(_datafile_fn_old)
        _datafiles_to_concatenate += _datafiles_fn_ranks

        if new_datafile_fn is None:
            _datafile_fn_new = os.path.join(_dir,'pyposmat.results.{}.out'.format(i_iteration))
        else:
            _datafile_fn_new = new_datafile_fn

        # now we need to merge the files
        _df = None

        for i,v in enumerate(_datafiles_to_concatenate):

            # an exception for the first iteration
            if i==0:
                datafile = PyposmatDataFile()
                datafile.read(filename=v)
                _df = datafile.df
                
            else:
                _iteration_id = '{}'.format(i_iteration)
                datafile = PyposmatDataFile()
                datafile.read(filename=v)
                _df = pd.concat([_df,datafile.df]).reset_index(drop=True)

        # create output file
        datafile = PyposmatDataFile()
        datafile.df = _df
        datafile.parameter_names = self.parameter_names
        datafile.qoi_names = self.qoi_names
        datafile.error_names = self.error_names

        if not os.path.exists(self.data_directory):
            os.mkdir(self.data_directory)
        
        try:
            datafile.write(filename=_datafile_fn_new)
        except FileNotFoundError as e:
            raise

 
    def merge_files(self,i_iteration,last_datafile_fn=None,new_datafile_fn=None):
        """ merge the pyposmat data files

        Args:
            i_iteration(int): the current iteration which just finished
            last_datafile_fn(str,optional): the filename of the last dataset in the data directory.
            new_datafile_fn(str,optional): where to output the file results 
        """

        if last_datafile_fn is None:
            last_datafile_fn = os.path.join(self.data_directory,
                                            'pyposmat.kde.{}.out'.format(i_iteration))

        if new_datafile_fn is None:
            new_datafile_fn = os.path.join(self.data_directory,
                                           'pyposmat.results.{}.out'.format(i_iteration))

        data_dir = self.data_directory
        rank_dirs = [v for v in os.listdir(self.root_directory) if v.startswith('rank_')]
        filenames = [os.path.join(self.root_directory,v,'pyposmat.results.out') for v in rank_dirs]

        data = None
        for i,v in enumerate(filenames):
            data_new = None
            if i == 0:
                data = PyposmatDataFile()
                data.read(filename=v)
            else:
                data_new = PyposmatDataFile()
                data_new.read(filename=v)

                data.df = pd.concat([data.df,data_new.df])

        nrows = len(data.df)
        sim_id_fmt = '{:0>2}_{:0>6}'
        sim_id_str = [sim_id_fmt.format(i_iteration,i) for i in range(nrows)]

        data.df['sim_id'] = [sim_id_fmt.format(i_iteration,i) for i in range(nrows)]

        data_old = PyposmatDataFile()
        try:
            data_old.read(filename=last_datafile_fn)
            data_old.df = pd.concat([data_old.df,data.df])
            data_old.write(filename=new_datafile_fn)
        except FileNotFoundError as e:
            if i_iteration == 0:
                data.write(filename=new_datafile_fn)
            else:
                raise


    def analyze_results(self,i_iteration,data_fn=None,config_fn=None,kde_fn=None):
        if data_fn is None:
            data_fn = os.path.join(\
                    self.root_directory,
                    self.data_directory,
                    'pyposmat.results.{}.out'.format(i_iteration))
        if config_fn is None:
            config_fn = os.path.join(\
                    self.root_directory,
                    self.configuration_filename)
        if kde_fn is None:
            kde_fn = os.path.join(\
                    self.root_directory,
                    self.data_directory,
                    'pyposmat.kde.{}.out'.format(i_iteration+1))

        data_analyzer = PyposmatDataAnalyzer()
        data_analyzer.read_configuration_file(filename=config_fn)
        data_analyzer.read_data_file(filename=data_fn)
        data_analyzer.write_kde_file(filename=kde_fn)

    @property
    def structure_directory(self):

        d = self.configuration.structures['structure_directory']

        if not os.path.isabs(d):
            d = os.path.join(self.root_directory,d)

        return d

    def read_configuration_file(self,filename=None):

        assert type(filename) in [type(None),str]
        assert type(self.configuration_filename) in [type(None),str]
        
        if filename is not None:
            self.configuration_filename = filename

        if not os.path.isabs(self.configuration_filename):
            self.configuration_filename = os.path.abspath(self.configuration_filename)

        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=self.configuration_filename)
        
        self.n_iterations = self.configuration.n_iterations
        self.qoi_names = self.configuration.qoi_names
        self.error_names = self.configuration.error_names
        self.parameter_names = self.configuration.parameter_names
  
        if self.mpi_rank == 0:
            self._write_parameter_names()
            self._write_qoi_names()
            self._write_error_names()

    def _write_parameter_names(self,parameter_names=None):
        if parameter_names is None: _parameter_names = self.parameter_names
        else: _parameter_names = parameter_names

        s = [80*'-']
        s += ['{:^80}'.format('PARAMETER_NAMES')]
        s += [80*'-']
        s += [p for p in _parameter_names]

        self.log("\n".join(s))
        
    def _write_qoi_names(self,qoi_names=None):
        if qoi_names is None: _qoi_names = self.qoi_names
        else: _qoi_names = qoi_names

        s = [80*'-']
        s += ['{:^80}'.format('QOI_NAMES')]
        s += [80*'-']
        s += [p for p in _qoi_names]

        self.log("\n".join(s))

    def _write_error_names(self,error_names = None):
        if error_names is None: _error_names = self.error_names
        else: _error_names = error_names

        s = [80*'-']
        s += ['{:^80}'.format('ERROR_NAMES')]
        s += [80*'-']
        s += [p for p in _error_names]

        self.log("\n".join(s))




if __name__ == "__main__":
    import Ni__eam__morse_exp_universal as Ni_eam

    #------------------------------------------------------------------------------
    # WRITE CONFIGURATION FILE
    #------------------------------------------------------------------------------
    Ni_eam_configuration = PyposmatConfigurationFile()
    Ni_eam_configuration.qois = Ni_eam.Ni_qoi_db.qois
    Ni_eam_configuration.potential = Ni_eam.Ni_eam_potential_formalism
    Ni_eam_configuration.structures = Ni_eam.Ni_structure_db
    Ni_eam_configuration.sampling_type = Ni_eam.Ni_eam_sampling
    Ni_eam_configuration.sampling_distribution =Ni_eam.Ni_eam_parameter_distribution
    Ni_eam_configuration.write(filename='pypospack.config.in')
    Ni_eam_configuration.read(filename='pypospack.config.in')

    pypospack_filename_in = 'pypospack.config.in'
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pypospack_filename_in)
    pyposmat_app.read_configuration_file()
    #pyposmat_app.read_configuration_file(filename=pyposmat_configuration_filename)
    pyposmat_app.run_all()
