from mpi4py import MPI 
import os, sys 

"""
https://gist.github.com/donkirkby/eec37a7cf25f89c6f259
"""

nproc = int(os.environ['SLURM_NTASKS'])
print("nprocs:{}".format(nproc))

mpi_info = MPI.Info.Create()
mpi_info.Set("add-hostfile", "worker_hosts")

comm = MPI.COMM_SELF.Spawn('vasp_std',
                           args=['>','vasp.log'],
                           maxprocs=nproc,
                           info=mpi_info).Merge()


