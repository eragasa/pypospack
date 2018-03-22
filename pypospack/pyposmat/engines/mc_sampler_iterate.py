import os,shutil,sys
import numpy as np
from mpi4py import MPI
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.engines import PyposmatMonteCarloSampler
from pypospack.pyposmat.data import PyposmatDataFile

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
                print(80*'-')
                print('{:80}'.format('BEGIN ITERATION {}/{}'.format(
                    i+1,self.n_iterations)))
                print(80*'-')
            MPI.COMM_WORLD.Barrier()

            self.run_simulations(i)
            MPI.COMM_WORLD.Barrier()

            if self.mpi_rank == 0:
                print('merging files')
                self.merge_files(i)
                self.analyze_results(i)
            MPI.COMM_WORLD.Barrier()
        print(80*'-')
        print('JOBCOMPLETE')

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

        # initialize()
        self.mc_sampler = PyposmatMonteCarloSampler(
                filename_in = _config_filename,
                filename_out = _results_filename,
                mpi_rank = self.mpi_rank,
                mpi_size = self.mpi_size)
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
            print(80*'-')
        MPI.COMM_WORLD.Barrier()

        _mc_config = self.mc_sampler.configuration.sampling_type[i_iteration]
        _mc_sample_type = _mc_config['type']
        _mc_n_samples = _mc_config['n_samples']

        # determine number of sims for this rank
        _n_samples_per_rank = int(_mc_n_samples/self.mpi_size)
        if _mc_n_samples%self.mpi_size > self.mpi_rank:
            _n_samples_per_rank += 1

        if _mc_sample_type == 'parametric':
            self.mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank)
        elif _mc_sample_type == 'kde':
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

        # get parameters from file
        elif _mc_sample_type == 'from_file':
            _filename_in = os.path.join(
                self.root_directory,
                _mc_config['file'])
            self.mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank,
                    filename=_filename_in
            )

        # return to root directory
        os.chdir(self.root_directory)

    def setup_mpi_environment(self):
        self.mpi_comm = MPI.COMM_WORLD
        self.mpi_rank = self.mpi_comm.Get_rank()
        self.mpi_size = self.mpi_comm.Get_size()
        self.mpi_procname = MPI.Get_processor_name()
        if self.mpi_rank == 0:
            self.print_mpi_environment()

    def print_mpi_environment(self):
        print(80*'-')
        print('{:^80}'.format('MPI COMMUNICATION INFORMATION'))
        print(80*'-')
        print('mpi_size={}'.format(self.mpi_size))

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
        
        print('Looking for previous kde file')
        print('    {}'.format(_filename_kde))

        
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
    
        print(self.parameter_names)
        print(self.qoi_names)
        print(self.error_names)

    def print_random_seeds(self):
            print(80*'-')
            print('{:^80}'.format('GENERATED RANDOM SEEDS'))
            print(80*'-')
            print()
            print('rv_seed={}'.format(self.rv_seed))
            print()
            print('{:^8} {:^8} {:^10}'.format('rank','iter','seed'))
            print('{} {} {}'.format(8*'-',8*'-',10*'-'))
            for i_rank in range(self.mpi_size):
                for i_iter in range(self.n_iterations):
                    print('{:^8} {:^8} {:>10}'.format(
                        i_rank,
                        i_iter,
                        self.rv_seeds[i_rank,i_iter]))
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
