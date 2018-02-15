import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatConfigurationFile

class PyposmatIterativeSampler(object):
    def __init__(self,
            configuration_filename):
        self.RANK_DIR_FORMAT = 'rank_{}'
        self.mpi_comm = None
        self.mpi_rank = None
        self.mpi_nprocs = None
        self.n_iterations = None
        self.rv_seed = None
        self.rv_seeds = None

        self.pyposmat_configuration_filename = configuration_filename
        self.pyposmat_configuration = None
        self.pyposmat_mc_sampler = None

        self.root_directory = os.getcwd()
        self.data_directory = 'data'

    def run_all(self):
        self.setup_mpi_environment()
        self.determine_rv_seeds()

        for i in range(self.n_iterations):
            if self.mpi_rank == 0:
                print(80*'-')
                print('{:80}'.format('BEGIN ITERATION {}/{}'.format(
                    i,self.n_iterations))
                print(80*'-')
       
            self.run_simulations(i)
            MPI.COMM_WORLD.Barrier()

            if self.mpi_rank == 0:
                self.merge_files(i)
                self.analyze_results(i)
            MPI.COMM_WORLD.Barrier()
    
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
                self.pyposmat_configuration_filename)

        _results_filename = os.path.join(
                self.root_directory,
                self.rank_directory,
                'pyposmat.results.out')

        # set random seed
        np.random.seed(self.rv_seeds[self.mpi_rank,i_iteration])
        
        self.pyposmat_mc_sampler = PyposmatMonteCarloSampler(
                filename_in = _config_filename,
                filename_out = _results_filename)
        self.pyposmat_mc_sampler.create_base_directories()
        self.pyposmat_mc_sampler.read_configuration_file()
        _structure_dir = self.pyposmat_mc_sampler.configuration.structures['structure_directory']
        self.pyposmat_mc_sampler.configuration.structures['structure_directory'] = \
                os.path.join('..',_structure_dir)
        self.pyposmat_mc_sampler.configure_qoi_manager()
        self.pyposmat_mc_sampler.configure_task_manager()
        self.pyposmat_mc_sampler.configure_pyposmat_datafile_out()
        #pyposmat_datafile_out = PyposmatDataFile(filename_out)

        if self.mpi_rank == 0:
            self.pyposmat_mc_sampler.print_structure_database()
            self.pyposmat_mc_sampler.print_sampling_configuration()
            self.pyposmat_mc_sampler.print_initial_parameter_distribution()

        _mc_config = self.pyposmat_mc_sampler.configuration.sampling_type[i_iteration]
        _mc_sample_type = _mc_config['type']
        _mc_n_samples = _mc_config['n_samples']
        
        # determine number of sims for this rank
        _n_samples_per_rank = int(_mc_n_samples/self.mpi_size)
        if _mc_n_samples%self.mpi_size > self.mpi_rank:
            _n_samples_per_rank += 1
        
        if _mc_sample_type == 'parametric':
            self.pyposmat_mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank)
        elif _mc_sample_type == 'kde':

            _filename_in = os.path.join(
                self.root_directory,
                self.data_directory,
                'pyposmat.kde.{}.out'.format(i_iteration))

            self.pyposmat_mc_sampler.run_simulations(
                    i_iteration=i_iteration,
                    n_samples=_n_samples_per_rank,
                    filename=_filename_in)
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
    
    def merge_files(self,i_iteration):
        n_ranks = self.mpi_size

        _data_directory = self.data_directory

        # delete directory if directory exists, then create it
        if i_iteration == 0:
            if os.path.isdir(self.data_directory):
                shutil.rmtree(self.data_directory)
            os.mkdir(self.data_directory)

        # gather all the lines into a list, then flush all at once
        str_list = []

        # grab the names, from the rank0 data file
        i_rank == 0:
        _rank_filename = os.path.join(
            'rank_{}'.format(i_rank),
            'pyposmat.results.out')
        with open(_rank_filename,'r') as f:
            lines = f.readlines()
        names = [v for v in lines[0].strip().split(',')]
        types = [v for v in lines[1].strip().split(',')]
        str_list.append(",".join(names))
        str_list.append(",".join(types))
        
        # first merge the kde file from the this iteration
        if i_iteration > 0:
            _filename_kde = os.path.join(
                self.data_directory,
                'pyposmat.kde.{}.out'.format(i_iteration))
            with open(_filename_kde,'r') as f:
                lines = f.readlines()
            for k,line in enumerate(lines):
                if k>1:
                    line = [v for v in lines[0].strip().split(',')]
                    str_list.append(",".join(line[len(names):]))
        
        # process all the results from this simulations
        for i_rank in range(n_ranks):
            # get the filename of the ith rank
            _rank_filename = os.path.join(
                'rank_{}'.format(i_rank),
                'pyposmat.results.out')
            print('...merging {}'.format(_rank_filename))

            with open(_rank_filename,'r') as f:
                lines = f.readlines()
            lines = [line.strip().split(',') for line in lines]
            for k,line in enumerate(lines):
                if k>1:
                    line[0] = '{}_{}_{}'.format(i_iteration,i_rank,line[0])
                    str_list.append(",".join(line))

        # write results to file
        _filename_out = os.path.join(
            self.data_directory,
            'pyposmat.results.{}.out'.format(i_iteration))
        with open(_filename_out,'w') as f_out:
            f_out.write("\n".join(str_list))
    
    def analyze_results(self,i_iteration):
        _filename_out = os.path.join(\
                self.root_directory,
                self.data_directory,
                'pyposmat.results.{}.out'.format(i_iteration)
        

    def read_configuration_file(self,filename=None):
        assert isinstance(filename,str) or filename is None

        _filename_in = None
        if filename is None:
            _filename_in = self.pyposmat_configuration_filename
        else:
            self.pyposmat_configuration_filename = filename
            _filename_in = filename
        
        self.pyposmat_configuration = PyposmatConfigurationFile(
                filename=_filename_in)
   
        self.n_iterations = self.pyposmat_configuration.sampling_type['n_iterations']

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
