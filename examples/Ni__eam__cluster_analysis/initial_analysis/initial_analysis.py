import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatClusterAnalysis
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import Pyposmat
pypospack_root = [v.strip() for v in os.environ['PYTHONPATH'].split(':') if v.endswith('pypospack')][0]

data_dir = os.path.join(pypospack_root,'data/Ni__eam__born_exp_fs__3.5NN')

configuration_fn = os.path.join(data_dir,'pyposmat.config.in')
datafile_fn = os.path.join(data_dir,'pyposmat.kde.0.out')

if __name__ == "__main__":

    d = OrderedDict([
        ('configuration_fn',configuration_fn),
        ('data_fn',datafile_fn),
        ('include_parameters',True),
        ('include_qois', False),
        ('include_errors',False)
    ])

    o = PyposmatClusterAnalysis.init_from_ordered_dict(d)

