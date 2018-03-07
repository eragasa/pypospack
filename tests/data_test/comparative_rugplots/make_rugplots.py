import copy
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
if __name__ == "__main__":
    _filename_pyposmat_config_fs  = 'subselect.d_metric.fs.out'
    _filename_pyposmat_config_bjs = 'subselect.d_metric.bjs.out'
    config=PyposmatConfigurationFile()
    config.read(filename='pyposmat.config.fs.in') 
     
    qoi_targets = OrderedDict()
    for k,v in config.qois.items():
        qoi_targets[k]=v['target']
    
    datafile_fs=PyposmatDataFile(filename=_filename_pyposmat_config_fs)
    datafile_fs.read()
    datafile_bjs=PyposmatDataFile(filename=_filename_pyposmat_config_bjs)
    datafile_bjs.read()

    error_names = datafile_fs.error_names
    qoi_names = datafile_fs.qoi_names
    
    for i_error,n_error in enumerate(error_names):
        _qoi_name = qoi_names[i_error]
        datafile_fs.df[n_error] = datafile_fs.df[n_error]/qoi_targets[_qoi_name]
        datafile_bjs.df[n_error] = datafile_bjs.df[n_error]/qoi_targets[_qoi_name]
    
    (_nrows_fs,_ncols) = datafile_fs.df.shape
    (_nrows_bjs,_ncols) = datafile_bjs.df.shape
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1,ncols=2)
    for i_error,n_error in enumerate(error_names):
        _yloc = [i_error+1]
        ax[0].scatter(
                datafile_fs.df[n_error],
                _nrows_fs*[i_error+1],
                marker='|',
                s=100.,
                color='k')
        ax[1].scatter(
                datafile_bjs.df[n_error],
                _nrows_fs*[i_error+1],
                marker='|',
                s=100.,
                color='k')
    plt.sca(ax[0])
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    plt.sca(ax[1])
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    fig.savefig('rugplots.png')

