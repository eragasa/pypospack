from pypospack.pyposmat.data import PyposmatDataFile
from collections import OrderedDict
filename_in ='pareto_009.out'
if __name__ == "__main__":
        _datafile_in = PyposmatDataFile(filename=filename_in)
        _datafile_in.read()

        mpi_size = 3
        mpi_rank = 0

        if mpi_rank is None:
            mpi_rank = 0
        else:
            mpi_rank = mpi_rank
        if mpi_size is None:
            mpi_size = 1
        else:
            mpi_size = mpi_size

        parameter_names = _datafile_in.parameter_names
        for row in _datafile_in.df.iterrows():
            if mpi_rank != row[0]%mpi_size:
                continue
            print(row[0])

                #_parameters = OrderedDict([(p,row[1][p]) for p in parameter_names])
                #print(_parameters)
