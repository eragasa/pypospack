from mpi4py import MPI

import numpy as np
import subprocess, sys, os, shutil

# Load the function to be evaluated
obj_dir = 'lmps_MgO_sample_single'
sys.path.insert(0, obj_dir)
from buckingham_single_as_a_func import obj_func


# Write some information to screen
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
sys.stdout.write('my rank is: '+ str(rank)+ ' and I live on: '+ name + '\n')
sys.stdout.write('the comm size is : '+ str(size) + '\n')


points = np.loadtxt('/home/jlarson/research/superchwm/data/input/MgO_buckingham/Random_sample_1e5')
points = points[:2,:]
num_pts = np.shape(points)[0]

# # ['chrg_Mg', 'chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_A', 'MgO_rho', 'MgO_C', 'OO_A', 'OO_rho', 'OO_C']
# low = np.array([  +1.5, 0, 0, 0.5,  0,  800, 0.29, 0,   500,  0.1, 25])      # Lower variable bounds
# upp = np.array([  +2.2, 0, 0, 0.5,  0, 1300, 0.33, 0, 25000,  0.4, 77])      # Upper variable bounds 
# np.random.seed(0)
# num_pts = 10000
# points = np.random.uniform(0,1,(num_pts,len(low)))*(upp-low)+low
# points[:,1] = -points[:,0] # The charge on O is the opposite of the charge on Mg



# Have each worker move to their own directory so function evaluations (e.g., LAMMPS files written) don't conflict.
worker_dir = '/scratch/' + obj_dir + "_" + str(rank)
# worker_dir = obj_dir + "_" + str(rank)
if os.path.exists(worker_dir):
	print("DELETING existing Worker directory.")
	sys.stdout.flush()
	shutil.rmtree(worker_dir)
saved_dir = os.getcwd()
shutil.copytree(obj_dir, worker_dir)
os.chdir(worker_dir)

# Perform the function evaluations
results = {}
for i in range(num_pts):
    if i % size == rank:
        # sys.stdout.write('rank ' + str(rank) + ' is working on point ' + str(i) + '\n')
        sys.stdout.write('rank ' + str(rank) + ' is working on point ' + str(i) + ' in directory' + os.getcwd() + '\n')
        results[i] = obj_func(points[i])


# Clean up
os.chdir(saved_dir)
shutil.rmtree(worker_dir)


# Save results to file
output_filename = 'LAMMPS_results_rank' + str(rank) + '.npy'
np.save(output_filename,results)

MPI.COMM_WORLD.Barrier()
# Have rank 0 clean everything up
if rank == 0:
    final_results = {}
    for i in range(0,size):
        rank_str = 'LAMMPS_results_rank' + str(i) + '.npy'
        tmp = np.load(rank_str).flatten()[0]
        final_results = {**final_results, **tmp}
        os.remove(rank_str)

    np.save('LAMMPS_results_5.npy',final_results)

