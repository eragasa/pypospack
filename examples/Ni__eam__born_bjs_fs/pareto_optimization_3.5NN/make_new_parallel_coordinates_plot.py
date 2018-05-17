import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


results_PunMishin2015 = OrderedDict([
    ('sim_id','PunMishin2015'),
    ('Ni_fcc.E_coh', -4.449999985713825),
    ('Ni_fcc.a0', 3.52000004514173),
    ('Ni_fcc.c11', 241.341629134211),
    ('Ni_fcc.c12', 150.824244634751),
    ('Ni_fcc.c44', 127.34413217099),
    ('Ni_fcc.B', 180.996706134571),
    ('Ni_fcc.G', 45.25869224972999),
    ('Ni_fcc.vac', 1.5722909444763218),
    ('Ni_fcc.100s', 0.12083783763315925),
    ('Ni_fcc.110s', 0.1306411055698218),
    ('Ni_fcc.111s', 0.1097944617790805),
    ('Ni_fcc.esf', 4.50280507300032e-14),
    ('Ni_fcc.isf', 0.008398360367557003),
    ('E_Ni_fcc_hcp', 0.02213171443965045),
    ('E_Ni_fcc_bcc', 0.06728584869762066),
    ('E_Ni_fcc_sc', 0.7235998449031951),
    ('E_Ni_fcc_dia', 1.4164289731208752)
])

results_Mishin1999 = OrderedDict([
    ('sim_id','Mishin1999'),
    ('Ni_fcc.E_coh', -4.449999998348875),
    ('Ni_fcc.a0', 3.51999943754043),
    ('Ni_fcc.c11', 247.862330912402),
    ('Ni_fcc.c12', 147.828379828392),
    ('Ni_fcc.c44', 124.838117591868),
    ('Ni_fcc.B', 181.17303018972868),
    ('Ni_fcc.G', 50.016975542005),
    ('Ni_fcc.vac', 1.601066770918635),
    ('Ni_fcc.100s', 0.11721226280222584),
    ('Ni_fcc.110s', 0.12813990050511612),
    ('Ni_fcc.111s', 0.10172841614230631),
    ('Ni_fcc.esf', -4.502806627505237e-14),
    ('Ni_fcc.isf', 0.007813859904327118),
    ('E_Ni_fcc_hcp', 0.020438536006599506),
    ('E_Ni_fcc_bcc', 0.11199380635944944),
    ('E_Ni_fcc_sc', 0.8330381769486048),
    ('E_Ni_fcc_dia', 0.010588663471387427)
])
results_Angelo1995 = OrderedDict([
    ('sim_id','Angelo1995'),
    ('Ni_fcc.E_coh', -4.4500000125938),
    ('Ni_fcc.a0', 3.52000035011041),
    ('Ni_fcc.c11', 246.715158886758),
    ('Ni_fcc.c12', 147.480621905903),
    ('Ni_fcc.c44', 124.960604806714),
    ('Ni_fcc.B', 180.55880089952134),
    ('Ni_fcc.G', 49.61726849042749),
    ('Ni_fcc.vac', 1.5938793445586157),
    ('Ni_fcc.100s', 0.1285743493399689),
    ('Ni_fcc.110s', 0.14781220704975925),
    ('Ni_fcc.111s', 0.12040311031473264),
    ('Ni_fcc.esf', 2.648708409963305e-15),
    ('Ni_fcc.isf', 0.005525903822777186),
    ('E_Ni_fcc_hcp', 0.012563864102899558),
    ('E_Ni_fcc_bcc', 0.07640412576407485),
    ('E_Ni_fcc_sc', 0.47501015461552987),
    ('E_Ni_fcc_dia', 0.012563864102899558)
])

ref_data = OrderedDict()
ref_data['PunMishin2015'] = results_PunMishin2015
ref_data['Mishin1999'] = results_Mishin1999
ref_data['Angelo1995'] = results_Angelo1995

ref_data_colors = OrderedDict()
ref_data_colors['PunMishin2015']= 'red'
ref_data_colors['Mishin1999'] = 'blue'
ref_data_colors['Angelo1995'] = 'green'

class PyposmatParallelCoordinates(object):
    pass

if __name__ == "__main__":
    # define the data directory
    data_directory = 'data'

    # read configuration file
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)

    # read the associated datafile
    datafile_fn = os.path.join(data_directory,'pyposmat.kde.5.out')
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_fn)

    plot_fn = 'parallelcoordinates_fs.png'

    #excluded_qoi_names = ['Ni_fcc.esf','Ni_fcc.isf','E_Ni_fcc_hcp']
    excluded_qoi_names = []
    qoi_names = [q for q in config.qoi_names if q not in excluded_qoi_names]
    print('qoi_names is length:{}'.format(len(qoi_names)))
    error_names = ["{}.err".format(q) for q in qoi_names]
    normed_error_names = ["{}.nerr".format(q) for q in qoi_names]
    qoi_targets = config.qoi_targets


    # calculate normalized error for sampled data
    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[nen] = datafile.df[qn]/q-1
    (nrows,ncols) = datafile.df.shape

    normederr_names = ['{}.nerr'.format(q) for q in qoi_names]
    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normederr_names]).sum(axis=1))
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from pandas.plotting import parallel_coordinates
    fig, ax = plt.subplots()


    reference_df = pd.DataFrame(list(ref_data.values()))
    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        reference_df[nen] = reference_df[qn]/q-1
    reference_df.set_index('sim_id')

    data_df = copy.deepcopy(datafile.df)
    data_df.set_index('sim_id')


    subselect_df = datafile.df.nsmallest(10,'d_metric')
    subselect_df.set_index('sim_id')

    is_plot_all_data = True
    if is_plot_all_data:
        parallel_coordinates(
                data_df[normed_error_names],
                'Ni_fcc.E_coh.nerr',
                color='grey'
                )

    parallel_coordinates(
            subselect_df[normed_error_names],
            'Ni_fcc.E_coh.nerr',
            color='k'
            )

    parallel_coordinates(
            reference_df[normed_error_names],
            'Ni_fcc.E_coh.nerr',
            )
    plt.gca().legend_.remove()

    #fig.savefig(plot_fn)
    plt.show()
    exit()
    (nr,nc)=df.shape
    print("We have {} potentials...".format(nr))
    for iqn,qn in enumerate(qoi_names):
        nen = '{}.nerr'.format(qn)
        x = df[nen]
        y = nr*[len(qoi_names)-iqn]
        ax.scatter(x,y,
            marker='|',
            #s=10.,
            color='grey')

    for ref_data_name,ref_data_dict in ref_data.items():
        for iqn,qn in enumerate(qoi_names):
            q = qoi_targets[qn]
            x = ref_data_dict[qn]/q-1
            y = len(qoi_names)-iqn
            ax.scatter(
                    x,
                    y,
                    s=300,
                    marker='|',
                    color=ref_data_colors[ref_data_name]
            )

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
    ax.legend(
        handles = [
            mpatches.Patch(
                color=ref_data_colors[ref_data_name],
                label=ref_data_name)
            for ref_data_name in ref_data
        ] + [
            mpatches.Patch(
                color='grey',
                label='best 50')
        ]
    )
    plt.show()
#    fig.savefig(plot_fn)
#    exit()
