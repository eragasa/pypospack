
class PyposmatIterativeSampler(object):
    def __init__(self,with_mpi=False):
        self.mpi_comm = None
        self.mpi_rank = None

        self.with_mpi = with_mpi

        if self.with_mpi is True:
            self.setup_mpi_environment()
        else:
            pass


    def setup_mpi_environment(self):
        from mpi4py import MPI

        self.mpi_comm = MPI.COMM_WORLD
        self.mpi_rank = self.mpi_comm.Get_rank()

    def set_sampling_seed(self):
        pass

if __name__ == "__main__":
    from mpi4py import MPI
    import numpy as np

    mpi_com = MPI.COMM_WORLD
    mpi_rank = mpi_com.Get_rank()
    mpi_nprocs = mpi_com.Get_size()
   
    print('hello, from rank {} of {} processors,'.format(mpi_rank,mpi_nprocs))
