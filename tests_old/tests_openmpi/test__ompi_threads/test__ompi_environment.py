import os
import sys
from mpi4py import MPI

def delete_old_files():
    files_to_remove = ['job.err','job.out','mpi_test.log']
    for fn in files_to_remove:
        if os.path.isfile('job.err'):
            os.remove(fn)

def get_mpi4py_version():
    pass

def test_thread_level():
    fn = "mpi_test.log"
    if MPI.COMM_WORLD.Get_rank() == 0:
        if not (MPI.Query_thread() == MPI.THREAD_MULTIPLE):
            with open(fn,'w') as f:
                m = "MPI does not provide enough thread support\n"
                m += 'MPI.Query_thread()={}\n'.format(MPI.Query_thread)
                f.write(m)
                print(m)
            sys.stderr.write("MPI does not provide enough thread support\n")
        else:
            with open(fn,'w') as f:
                m = 'MPI.Query_thread()={}\n'.format(MPI.Query_thread)
                f.write(m)
                print(m)
    #sys.exit(0)

if __name__ == "__main__":
    delete_old_files()
    test_thread_level()
    #comm = MPI.COMM_WORLD
    #mpi_rank = comm.Get_rank()
    #if mpi_rank == 0:
    #    with open(fn,'w') as f:
    #        print(MPI.Query_thread())
    #        f.write(str(MPI.Query_thread()))
