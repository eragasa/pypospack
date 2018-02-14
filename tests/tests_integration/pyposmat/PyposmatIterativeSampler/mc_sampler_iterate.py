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
        self.pyposmat_mc_sampler = None

        self.root_directory = os.getcwd()
        
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

        # determine number of sims for this rank
        _n_sims = 1000 * self.mpi_size +2
        _n_sims_per_rank = int(_n_sims/self.mpi_size)
        if _n_sims%self.mpi_size > self.mpi_rank:
            _n_sims_per_rank += 1
        
        # set random seed
        np.random.seed(self.rv_seeds[self.mpi_rank,i_iteration])
        
        self.pyposmat_mc_sampler = PyposmatMonteCarloSampler(
                filename_in = _config_filename,
                filename_out = _results_filename)
        #self.pyposmat_mc_sampler.create_base_directories()
        #self.pyposmat_mc_sampler.read_configuration_file()
        #self.pyposmat_mc_sampler.configure_qoi_manager()
        #self.pyposmat_mc_sampler.pyposmat_datafile_out = PyposmatDataFile(filename_out)

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
    
    def run_all(self):
        self.setup_mpi_environment()
        self.determine_rv_seeds()

        for i in range(self.n_iterations):
            self.run_simulations(i)
            MPI.COMM_WORLD.Barrier()

            self.analyze_results(i)
            MPI.COMM_WORLD.Barrier()
    def analyze_results(self,i_iterations):
        # merge datafiles
        if self.mpi_rank == 0:
            _filename_out = 'pyposmat.results.{}.out'.format(
                    i_iterations)
            print('merging results into {}'.format(_filename_out))
            for i_rank in range(self.mpi_size):
                _filename_in = os.path.join(
                    self.RANK_DIR_FORMAT.format(i_rank))
                print('\t merging...{}'.format(_filename_in))

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

    pypospack_filename_in = 'pypospack.config.in'
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pypospack_filename_in)
    #pyposmat_app.read_configuration_file(filename=pyposmat_configuration_filename)
    pyposmat_app.n_iterations = 10
    pyposmat_app.run_all()
