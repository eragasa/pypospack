import os
data_dir = "Ni__eam__morse_exp_bjs_00"

def analyze_data_directory(data_dir):
    i = 0
    contents = []
    if not os.path.exists(data_dir): return i, contents
    if not os.path.isdir(data_dir): return i, contents

    while True:
        kde_fn = os.path.join(data_dir,"pyposmat.kde.{}.out".format(i))
        if os.path.exists(kde_fn):
            contents.append(kde_fn)
        else:
            if i > 0:
                contents.append(results_fn)
                break

        results_fn = os.path.join(data_dir,"pyposmat.results.{}.out".format(i))
        if os.path.exists(results_fn): pass
        else:break

        i = i + 1

    return i, contents

root_dir = "Ni__eam__born_exp_bjs_00" 

def analyze_rank_directories(root_dir):
    i = 0
    contents = []

    while True:
        rank_dir = os.path.join(root_dir,"rank_{}".format(i))
        if not os.path.exists(rank_dir): break
        if not os.path.isdir(rank_dir): break
        rank_fn = os.path.join("rank_{}".format(i),"pyposmat.results.out")
        if not os.path.exists(os.path.join(root_dir,rank_fn)): break
        if not os.path.isfile(os.path.join(root_dir,rank_fn)):break
        else: contents.append(rank_fn)
        i = i + 1
    return i,contents

from pypospack.pyposmat.data import PyposmatDataFile
import pandas as pd
def merge_pypospack_datafiles(datafile_fns):
    d0 = PyposmatDataFile()
    d0.read(filename=datafile_fns[0])
    df0 = d0.df
    for i in range(1,len(datafile_fns)):
        print("merging {}...".format(datafile_fns[i]))
        d = PyposmatDataFile()
        d.read(filename=datafile_fns[i])
        df = d.df
        
        df0 = pd.concat([df0,df]).drop_duplicates().reset_index(drop=True)
    d0.df = df0
    return d0
    
i_iterations, data_contents = analyze_data_directory(data_dir)
n_ranks, rank_contents = analyze_rank_directories(root_dir)
df = merge_pypospack_datafiles(data_contents)
df.write("pypospack.recovery.out")
contents = data_contents + rank_contents

print("i_iterations:{}".format(i_iterations))
print("n_ranks:{}".format(n_ranks))
