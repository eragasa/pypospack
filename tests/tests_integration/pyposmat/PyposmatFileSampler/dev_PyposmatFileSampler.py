import pytest
import os, copy
from collections import OrderedDict
import numpy as np
import pandas as pd
from pypospack.pyposmat import PyposmatDataFile

data_directory = os.path.join(
        '../../../test_data/test_PyposmatData',
        'data',
        'output')
pypospack_data_filename = os.path.join(
        data_directory,
        'culled_009.out')
names = [
        'chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C',
        'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
        'MgO_NaCl.B','MgO_NaCl.G',
        'MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch',
        'MgO_NaCl.001s',
        'MgO_NaCl.a0.err','MgO_NaCl.c11.err','MgO_NaCl.c12.err',
        'MgO_NaCl.c44.err','MgO_NaCl.B.err','MgO_NaCl.G.err',
        'MgO_NaCl.fr_a.err','MgO_NaCl.fr_c.err','MgO_NaCl.sch.err',
        'MgO_NaCl.001s.err']
latex_labels = {
        'MgO_NaCl.a0':r'$a_0$',
        'MgO_NaCl.c11':r'$c_{11}$',
        'MgO_NaCl.c12':r'$c_{12}$',
        'MgO_NaCl.c44':r'$c_{44}$',
        'MgO_NaCl.B':r'$B$',
        'MgO_NaCl.G':r'$G$',
        'MgO_NaCl.fr_a':r'$E_{fr,a}$',
        'MgO_NaCl.fr_c':r'$E_{fr,c}$',
        'MgO_NaCl.sch':r'$E_{sch}$',
        'MgO_NaCl.001s':r'$\gamma_{001}$'}
parameter_names = ['chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C']
qoi_names = [
        'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
        'MgO_NaCl.B','MgO_NaCl.G', 'MgO_NaCl.fr_a','MgO_NaCl.fr_c',
        'MgO_NaCl.sch','MgO_NaCl.001s']
error_names = ['{}.err'.format(s) for s in qoi_names]
qoi_reference_dft = {
        'MgO_NaCl.a0': 4.246,
        'MgO_NaCl.c11': 277.00031,
        'MgO_NaCl.c12': 91.67016,
        'MgO_NaCl.c44': 144.00722,
        'MgO_NaCl.B': 153.4468767,
        'MgO_NaCl.G': 92.665075,
        'MgO_NaCl.fr_a': 10.9781666,
        'MgO_NaCl.fr_c': 8.98642095,
        'MgO_NaCl.sch':5.067179685,
        'MgO_NaCl.001s': 0.055950069}
reference_LC = [2.0,-2.0,0.0,0.5,0.0,821.6,0.3242,0.0,22764.0,0.149,27.88,4.21078276561128,307.5718095128,171.135602331774,168.168424521017,216.61433805878266,68.21810359051298,9.68024989411606,9.810715180656189,5.797051474639375,0.06783817649966246,-0.035217234388720264,30.571809512799973,79.46560233177401,24.158424521016997,63.164338058782675,-24.441896409487015,-1.2977501058839405,0.8247151806561881,0.7300514746393745,0.011888176499662464]
reference_BG1 = [2.0,-2.0,0.0,0.5,0.0,1279.69,0.29969,0.0,9547.96,0.21916,32.0,4.20923604431415,383.274119165401,169.434215310753,179.601185701851,240.71418326230233,106.91995192732399,12.419259511088967,11.869114175328832,7.198887069605007,0.08070791160146304,-0.036763955685850114,106.27411916540098,77.76421531075299,35.591185701851,87.26418326230234,14.259951927323996,1.4412595110889672,2.8831141753288314,2.131887069605007,0.02475791160146304]
reference_BG2 = [1.7,-1.7,0.0,0.5,0.0,929.69,0.29909,0.0,4870,0.2679,77.0,4.222448,301.315822490901,150.827961179874,142.471471673523,200.990581616883,75.2439306555135,10.434727086962994,8.526633932683126,5.509135247188169,0.0692527749868838,-0.02355200000000046,24.315822490900985,59.15796117987399,-1.5385283264769782,47.540581616883,-17.416069344486502,-0.543272913037006,-0.4593660673168749,0.442135247188169,0.0133027749868838]


class PyposmatFileSampler(object):

    def __init__(self,
            pyposmat_data_filename):
        self.pyposmat_data_filename = pyposmat_data_filename
        self.workflow = None
        self.tasks = None
        self.obj_pyposmat_datafile = None

    def add_task(self,task_name,task):
        if self.tasks is None:
            self.tasks = OrderedDict()
        self.task[task_name] = task

    def read_pyposmat_datafile(self,filename=None):
        if filename is not None:
            self.pyposmat_data_filename= filename

        self.obj_pyposmat_datafile = PyposmatDataFile(
                filename=self.pyposmat_data_filename)

    def evaluate_from_file(self,filename=None):
        if filename is not None:
            self.pypospack_data_filename = filename
            self.read_pyposmat_datafile()
        if self.obj_pyposmat_datafile is None:
            self.read_pyposmat_datafile():

        _dataframe = self.obj_pypospack_datafile.df
        self.obj_pypospack_datafile = PyposmatDataFile(
                filename = self.filename)
        self.obj_pypospack_datafile.read(
                filename = self.filename)
        for row in self.obj_pypospack_datafile.df:
            self.evaluate_parameter_set(self,parameters)

    def evaluate_parameter_set(self,parameters):
        print(parameters)
if __name__ == "__main__":
    n = 1000
    scaling_factors='DFT'
    err_type='abs'
    pypospack_optimal_filename='optimal_{n}.out'.format(n=n)

    obj_datafile = PyposmatDataFile(filename=pypospack_data_filename)
    obj_datafile.read(filename=pypospack_data_filename)

    obj_datafile.qoi_references = OrderedDict()
    obj_datafile.qoi_references['DFT'] = copy.deepcopy(qoi_reference_dft)
    obj_datafile.create_optimal_population(
            n=n,
            scaling_factors=scaling_factors,
            err_type=err_type)

    obj_datafile.write_optimal_population(
            filename=pypospack_optimal_filename,
            n=n,
            scaling_factors=scaling_factors,
            err_type=err_type)
    obj_datafile.read(filename=pypospack_optimal_filename)
