import copy,os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":
    data_directory = os.path.join('../../../../',
            'data_test',
            'Ni__eam__born_exp_fs_00',
            'data__Ni__eam__born_exp_fs_02')
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)

    datafile_fn = os.path.join(data_directory,'pyposmat.kde.9.out')
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_fn)

    error_names = config.error_names
    qoi_names = config.qoi_names
    qoi_targets = config.qoi_targets

    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[nen] = datafile.df[qn]/q-1
   
    (nrows,ncols) = datafile.df.shape
    
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    (nrows,ncols) = datafile.df.shape
    for iqn,qn in enumerate(qoi_names):
        nen = '{}.nerr'.format(qn)
        x = datafile.df[nen]
        y = nrows*[len(qoi_names)-iqn]
        ax.scatter(x,y,
            marker='|',
            s=100.,
            color='k')
 
    ax.set_ylim([0,len(qoi_names)+1])
    ax.set_xlabel('Pct Error Difference')
    yticks_labels = [config.latex_labels[qn]['name'] for qn in qoi_names]
    plt.sca(ax)
    plt.yticks(
            list(range(1,len(qoi_names)+1)),
            list(reversed(yticks_labels))
        )
    fig.savefig('rugplots.eps')
    
    exit()
