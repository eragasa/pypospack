import copy,yaml
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatMonteCarloSampler
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
filename_out='pyposmat.results.out'
pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
print('pypospack_root_directory:{}'.format(pypospack_root_dir))
#--------------------------------------------------------------------------------------------
# WRITE CONFIGURATION FILE
#--------------------------------------------------------------------------------------------
data_directory=os.path.join(
        pypospack_root_dir,
        'data/Ni__eam__born_exp_rose/preconditioning_3.5NN'
Ni_eam_configuration = PyposmatConfigurationFile()

