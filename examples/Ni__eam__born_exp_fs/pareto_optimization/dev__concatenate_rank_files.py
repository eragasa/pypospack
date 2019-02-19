from pypospack.pyposmat.data import PyposmatDataFile
import os
import pandas as pd

i_iteration = 1
data_dir = "./"
rank_dir = [v for v in os.listdir(data_dir) if v.startswith('rank_')]
filenames = [os.path.join(data_dir,v,'pyposmat.results.out') for v in rank_dir]
data_old_fn = 'pyposmat.kde.old'
data = None
for i,v in enumerate(filenames):
    data_new = None
    if i == 0:
        data = PyposmatDataFile()
        data.read(filename=v)
    else:
        data_new = PyposmatDataFile()
        data_new.read(filename=v)

        data.df = pd.concat([data.df,data_new.df])

nrows = len(data.df)
sim_id_fmt = '{:0>2}_{:0>6}'
sim_id_str = [sim_id_fmt.format(i_iteration,i) for i in range(nrows)]

data.df['sim_id'] = [sim_id_fmt.format(i_iteration,i) for i in range(nrows)]
print(sim_id_str)
print(data.df)

data_old = PyposmatDataFile()
data_old.read(filename=data_old_fn)

data_old.df = pd.concat([data_old.df,data.df])
data_old.write(filename='pyposmat.data.new')
