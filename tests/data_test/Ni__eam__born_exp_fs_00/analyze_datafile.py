import copy
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
if __name__ == "__main__":
    _filename_pyposmat_config = 'pyposmat.config.in'
    config=PyposmatConfigurationFile()
    config.read(filename='pyposmat.config.in') 
     
    qoi_targets = OrderedDict()
    for k,v in config.qois.items():
        qoi_targets[k]=v['target']
    
    _filename_pyposmat_data = 'data/pyposmat.kde.10.out'
    _n_potentials = 1000
    datafile=PyposmatDataFile(filename=_filename_pyposmat_data)
    datafile.read()
    datafile.qoi_references = OrderedDict()
    datafile.qoi_references['TARGET'] = copy.deepcopy(qoi_targets)
    datafile.score_by_d_metric(scaling_factors='TARGET')
    datafile.subselect_by_score(
            score_name='d_metric',
            n=_n_potentials)
    _filename_subselect = datafile.write_subselect()

    datafile=PyposmatDataFile(filename=_filename_subselect)
    datafile.read()
    error_names = datafile.error_names
    qoi_names = datafile.qoi_names
    for i_error,n_error in enumerate(error_names):
        _qoi_name = qoi_names[i_error]
        datafile.df[n_error] = datafile.df[n_error]/qoi_targets[_qoi_name]

    (_nrows,_ncols) = datafile.df.shape
    import matplotlib.pyplot as plt
    print(datafile.df[error_names])
    fig, ax = plt.subplots()
    for i_error,n_error in enumerate(error_names):
        _yloc = [i_error+1]
        ax.scatter(
                datafile.df[n_error],
                _nrows*[i_error+1],
                marker='|',
                s=100.,
                color='k')
    plt.sca(ax)
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    fig.savefig('rugplots.png')

