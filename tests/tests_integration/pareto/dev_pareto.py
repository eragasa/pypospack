import os
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
#import pareto

from random     import randint, seed
import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('%s function took %0.3f ms' % (f.__name__, (time2-time1)*1000.0))
        return ret
    return wrap

def dominates(p1, p2):
    for x,y in zip(p1,p2):
        if y < x:
            return False
    return True

def pareto_bruteforce(pts, indices = None):
    if indices is None:
        indices = list(range(len(pts)))
    result = []
    for i in indices:
        for j in indices:
            if i == j: continue
            if dominates(pts[j], pts[i]):
                break
        else:
            result.append(i)
    return result

def pareto_merge(lo, hi, i, dim):
    if len(lo) == 0 or len(hi) == 0:
        return lo + hi

    survivors = set()
    for j in range(dim):
        if i == j: continue
        m = min(p[j] for p in lo)
        survivors.update(k for k in range(len(hi)) if hi[k][j] < m)

    return lo + [hi[k] for k in survivors]

def pareto(pts, indices = None, i = 0):
    if indices is None:
        indices = list(range(len(pts)))
    l = len(indices)
    if l <= 1:
        return indices

    if l < 1000:
        return pareto_bruteforce(pts, indices)

    indices.sort(key = lambda x: pts[x][i])     # lazy: should use partition instead

    dim = len(pts[0])
    optimalLo = pareto(pts, indices[:l//2], (i + 1) % dim)
    optimalHi = pareto(pts, indices[l//2:], (i + 1) % dim)

    return pareto_bruteforce(pts, optimalLo + optimalHi)     # lazy: FIXME
    #return pareto_merge(optimalLo, optimalHi, i, dim)

def read_data(fn):
    with open(fn) as f:
        lines = f.readlines()
        names = ['sim_id'] + [x.split() for x in lines[0].split('|')]
        names[1] = names[1][1:]     # skip sim_id
        values = []
        for line in lines[1:]:
            line = line.split('|')
            line = [x.split() for x in line]
            values.append([
                int(line[0][0]), 
                [float(x) for x in line[0][1:]], 
                [float(x) for x in line[1]]])

        return names, values

class PypospackDataAnalyzer(object):
    def __init__(self):
        self._pyposmat_configuration = None

    @property
    def pyposmat_configuration(self):
        return self._pypospack_configuration

    @pyposmat_configuration.setter
    def pypospack_configuration(self,configuration):
        assert type(configuration) is PypospackConfigurationFile

    @property
    def parameter_names(self): 
        return self._pyposmat_datafile.parameter_names

    @property
    def qoi_names(self):
        return self._pyposmat_datafile.qoi_names

    @property
    def error_names(self):
        return self._pyposmat_datafile.error_names

    def read_pypospack_configuration_file(self,filename):
        self._pyposmat_configuration = PyposmatConfigurationFile()
        self._pyposmat_configuration.read(filename)

    def read_pypospack_datafile(self,filename):
        self._pyposmat_datafile(filename)
        self._pyposmat_datafile.read()

        self._df = copy.deepcopy(self._pyposmat_datafile.df)
        self._df['sim_id'] = df['sim_id'].apply(np.int64)
        self._df[self.error_names] = self._df[self.error_names].abs()

    def calculate_pareto_set(self):
        _results = []
        for i_row, row in self._df.iterrows():
            _results.append([
                    int(row['sim_id']),
                    row[self.parameter_names].values.tolist(),
                    row[self.error_names].values.tolist()
                    ])

        is_pareto_idx = pareto(_results)    

if __name__ == "__main__":
    data_directory = 'data'
    pyposmat_data_filename = 'pypospack.results.out'
   
    datafile = PyposmatDataFile(filename=os.path.join(
        data_directory,pyposmat_data_filename))

    datafile.read()

    print(datafile.parameter_names)
    print(datafile.qoi_names)
    print(datafile.error_names)
    #print(datafile.df)

    import copy
    import numpy as np
    df = copy.deepcopy(datafile.df)

    p_names = datafile.parameter_names
    q_names = datafile.qoi_names
    e_names = datafile.error_names

    p_column_idx = [df.columns.get_loc(n) for n in p_names if n in df.columns]
    q_column_idx = [df.columns.get_loc(n) for n in q_names if n in df.columns]
    e_column_idx = [df.columns.get_loc(n) for n in e_names if n in df.columns]

    df['sim_id'] = df['sim_id'].apply(np.int64)
    df[e_names] = df[e_names].abs()

    values = []
    i_iter = 0
    for i_row, row in df.iterrows():
        values.append(
                [
                    "{}_{}".format(i_iter,int(row['sim_id'])),
                    row[p_column_idx].values.tolist(),
                    row[e_column_idx].values.tolist()
                ])
    for idx,p,v in values:
        print('idx:',idx)
        print('p:',p)
        print('v:',v)
    
    is_pareto_idx = pareto([v for idx,p,v in values])
    df = copy.deepcopy(datafile.df)
    df['is_pareto'] = 0
    df.loc[is_pareto_idx,'is_pareto'] = 1
    n_results = len(values)
    n_pareto = len(is_pareto_idx)
    print('n_results={}'.format(n_results))
    print('n_pareto={}'.format(n_pareto))

    kde_parameters = df[df['is_pareto'] == 1][p_names]
    print(kde_parameters)
