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
import os
from pypospack.visualization import ParetoOptimizationVisualization


import numpy as np
import pandas as pd
class PypospackVisualization(ParetoOptimizationVisualization):

    def load_data_file(self, fname):
        names = []
        names_types = []
        data = []

        # read in the line
        lines = None
        try:
            with open(fname, 'r') as f:
                lines = f.readlines()
        except:
            raise

        for i, line in enumerate(lines):
            line = line.strip()
            if i == 0:
                names = [v.strip() for v in line.split(',')]
            elif i == 1:
                name_types = [v.strip() for v in line.split(',')]
            else:
                data_line = [v.strip() for v in line.split(',')]
                for j in range(1,len(data_line)):
                    data_line[j] = float(data_line[j])
                data.append(data_line)
       
        #data = np.array(data)

        assert len(names) == len(name_types)
        # block below organizes data names by type into lists
        param_names = []
        param_key_index = []
        qoi_names = []
        qoi_key_index = []
        err_names = []
        err_key_index = []
        for i, v in enumerate(name_types):
            if v == 'param':
                param_names.append(names[i])
                param_key_index.append(i)
            elif v == 'qoi':
                qoi_names.append(names[i])
                qoi_key_index.append(i)
            elif v == 'err':
                err_names.append(names[i])
                err_key_index.append(i)

        
        self.df = pd.DataFrame(data,columns=names)
       
        _names = [n for n in names if n is not 'sim_id']
        self.total_df = self.df[_names]
        self.param_df = self.df[param_names]
        self.qoi_df = self.df[qoi_names]
        self.err_df = self.df[err_names]
        
        self.param_names = list(param_names)
        self.qoi_names = list(qoi_names)
        self.err_names = list(err_names)

if __name__ == "__main__":
    data_dir = 'data__morse_exp_fs_2'
    filename = os.path.join(data_dir,'pyposmat.results.9.out')

    vizdemo = PypospackVisualization()
    vizdemo.load_data_file(fname= filename)
    print(80*'-')
    print('{:^80}'.format('PARAMETERS'))
    print(80*'-')
    for v in vizdemo.param_names:
        print(v)
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
