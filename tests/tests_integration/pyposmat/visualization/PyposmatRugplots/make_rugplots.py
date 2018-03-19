import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":
    data_directory = os.path.join('../../../../',
            'data_test',
            'Ni__eam__born_exp_bjs_00',
            'data__Ni__eam__born_exp_bjs_04i')
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)

    datafile_fn = os.path.join(data_directory,'pyposmat.kde.20.out')
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_fn)

    qoi_names = [q for q in config.qoi_names if q not in ['Ni_fcc.esf','Ni_fcc.isf','E_Ni_fcc_hcp']]
    print('qoi_names is length:{}'.format(len(qoi_names)))
    error_names = ["{}.err".format(q) for q in qoi_names]
    qoi_targets = config.qoi_targets

    
    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[nen] = datafile.df[qn]/q-1
   
    (nrows,ncols) = datafile.df.shape
    
    normederr_names = ['{}.nerr'.format(q) for q in qoi_names]
    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normederr_names]).sum(axis=1))
    df = datafile.df.nsmallest(50,'d_metric')
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    (nr,nc)=df.shape
    print("We have {} potentials...".format(nr))
    for iqn,qn in enumerate(qoi_names):
        nen = '{}.nerr'.format(qn)
        x = df[nen]
        y = nr*[len(qoi_names)-iqn]
        ax.scatter(x,y,
            marker='|',
            #s=10.,
            color='k')
    plt.axvline(0,color='k',linestyle='-',linewidth=.1)
    ax.set_xlabel('Pct Error Difference')
    yticks_labels = [config.latex_labels[qn]['name'] for qn in qoi_names]
    print("length of yticks_labels:{}".format(len(yticks_labels)))
    plt.sca(ax)
    plt.yticks(
            list(range(1,len(qoi_names)+1)),
            list(reversed(yticks_labels))
        )
    #vals = ax.get_xticks()
    #ax.set_xticklabels(
    #        ['{:.1%}'.format(int(x*100)) for x in vals]
    #    )
    ax.set_xlim([-1,1])
    ax.set_ylim([0,len(qoi_names)+1])
    fig.savefig('rugplots_bjs.eps')
    
    exit()
