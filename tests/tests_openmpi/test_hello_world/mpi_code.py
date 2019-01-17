from mpi4py import MPI

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    with open("rank_{}".format(rank), "w") as f:
        f.write("success")
