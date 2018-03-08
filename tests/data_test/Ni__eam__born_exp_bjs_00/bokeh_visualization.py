# imports here
'''
TODO:
scaling and box zoom

probability density plot
-histogram
1d fit
2s fit
-multivariant normal distribution

kernel density estimate
-pyflames post module

principal components analysis
clustering
manifold learning T-sne
'''
import os,copy
from pypospack.visualization import ParetoOptimizationVisualization
from pypospack.pyposmat.data import PyposmatDataFile

import numpy as np
import pandas as pd
class PypospackVisualization(ParetoOptimizationVisualization):

    def load_data_file(self, fname):
        data = PyposmatDataFile().read(fname)
        self.param_names = list(data.parameter_names)
        self.qoi_names = list(data.qoi_names)
        self.err_names = list(data.err_names)
        self.names = list(
                data.parameter_names
                + data.qoi_names
                + data.err_names)
        self.df = copy.deepcopy(data.df)
        self.total_df = self.df[self.names]
        self.param_df = self.df[self.param_names]
        self.qoi_df = self.df[self.qoi_names]
        self.err_df = self.df[self.err_names]

        self.param_names = list(param_names)
        self.qoi_names = list(qoi_names)
        self.err_names = list(err_names)

if __name__ == "__main__":
    #data_dir = 'data'
    #filename = os.path.join(data_dir,'pyposmat.kde.5.out')
    filename='data__Ni__eam__born_exp_bjs_00/pyposmat.kde.10.out'
    vizdemo = PypospackVisualization()
    vizdemo.load_data_file(fname= filename)
    print(80*'-')
    print('{:^80}'.format('PARAMETERS'))
    print(80*'-')
    for v in vizdemo.param_names:
        print(v)
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
