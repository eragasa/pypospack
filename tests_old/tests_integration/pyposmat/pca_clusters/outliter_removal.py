import os
import numpy as np
import pandas as pd
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from scipy.stats import gaussian_kde
from pypospack.kde import Chiu1999_h

def filter_by_znormalized_errors(df,percentile,qoi_names):
    for qn in qoi_names:
        en = "{}.err".format(qn)
        zen = "{}.zerr".format(qn)
        df[zen] = (df[en]-df[en].mean())/df[en].std()

    zerror_names = ["{}.zerr".format(q) for q in qoi_names]
    df['z_error_dist'] = np.sqrt(np.square(df[zerror_names]).sum(axis=1))

    nr0,nc0 = df.shape
    nr1 = int(percentile*nr0//1)
    df=df.nsmallest(nr1,"z_error_dist").reset_index(drop=True)

    return df

if __name__ == "__main__":
    resource_dir = 'resources'
    config_fn = os.path.join(resource_dir,'pyposmat.config.in')
    data_fn = os.path.join(resource_dir,'pyposmat.kde.14.out')

    percentile = .90
    config = PyposmatConfigurationFile()
    config.read(filename=config_fn)

    data = PyposmatDataFile()
    data.read(filename=data_fn)

    qoi_names = config.qoi_names

    df = filter_by_znormalized_errors(data.df,percentile,qoi_names)
    _df = data.df
    for qn in qoi_names:
        en = "{}.err".format(qn)
        z_en = "{}.zerr".format(qn)
        _df[z_en] = (_df[en])/_df[en].std()

    zerror_names = ["{}.zerr".format(q) for q in qoi_names]
    _df['z_err_dist'] = np.sqrt(np.square(_df[zerror_names]).sum(axis=1))
    nr0,nc0 = _df.shape

    print("total_number_of_rows:{}".format(nr0))
    nr1 = int(percentile*nr0//1)
    print("keeping_number_of_rows:{}".format(nr1))
    _df=_df.nsmallest(nr1,'z_err_dist').reset_index(drop=True)
    nr2,nc2 = _df.shape
    print("current_df_size:{},{}".format(nr2,nc2))
    data.df = _df.copy(deep=True)
    nr3,nc3 = data.df.shape
    print("current_df_size:{},{}".format(nr3,nc3))
    import matplotlib.pyplot as plt
    plt.figure()
    _df['z_err_dist'].plot.hist()
    plt.show()
