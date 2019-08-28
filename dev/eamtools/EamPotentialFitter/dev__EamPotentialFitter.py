import os
import pypospack.utils
from potentialfitter import EamPotentialFitter

def initialize_from_setfl_file(filename):
    o = SeatonSetflReader(path=filename)
    print(symbols)

if __name__ == "__main__":
    filename = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','potentials','Ni__eam',
            ' Mishin-Ni-Al-2009.eam.alloy') 

    initialize_from_setfl_file(filename=filename)




