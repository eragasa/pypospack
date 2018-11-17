import os,sys
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

if __name__ == "__main__":
    src_dir='data'
    dst_dir='analysis_Si_sw_vac'
    
    if not os.path.exists(dst_dir):
        os.mkdir(dst_dir)

    config_fn=os.path.join(src_dir,'pyposmat.config.in')
    print("reading the configuration file {}...".format(
        config_fn))
    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)

    n_iterations = config.n_iterations
    print("there are {} iterations".format(
        n_iterations))

    # LOAD DATAFILES ----------------------------------------------------------
    results_data=[]
    kde_data=[]
    for i in range(n_iterations):
        kde_fn = os.path.join(
                src_dir,
                "pyposmat.kde.{}.out".format(i))
        results_fn = os.path.join(
                src_dir,
                "pyposmat.results.{}.out".format(i))
       
        print("reading {}...".format(kde_fn))
        try:
            kde_data.append(PyposmatDataFile())
            kde_data[i].read(filename=kde_fn)
        except FileNotFoundError as e:
            kde_data[i] = None

        #print("reading {}...".format(results_fn))
        #try:
        #    results_data.append(PyposmatDataFile())
        #    results_data[i].read(filename=results_fn)
        #except FileNotFoundError as e:
        #    results_data[i] = None

    N_pareto_pts=[]
    dN_pareto_pts=[]
    for i in range(n_iterations):
        if i == 0:
            nr0=0
            nc0=0
            if kde_data[i] is not None:
                nr1,nc1 = kde_data[i].df.shape
            else:
                nr1=0
                nc1=0
        else:
            if kde_data[i-1] is None:
                nr0=0
                nc0=0
            else:
                nr0,nc0 = kde_data[i-1].df.shape
            nr1,nc1 = kde_data[i].df.shape

        dN_pareto_pts.append(nr1-nr0)
        N_pareto_pts.append(nr1)

    for i in range(n_iterations):
        print("{:10} {:10}".format(
            dN_pareto_pts[i],
            N_pareto_pts[i]))

    for q in config.qoi_names:
        aen = "{}.abserr".format(q)
        en = "{}.err".format(q)
        for i in range(n_iterations):
            kde_data[i].df[aen] = kde_data[i].df[en].abs()
        config.latex_labels[aen] = OrderedDict()
        config.latex_labels[aen]['name'] = r"$|$"+config.latex_labels[q]['name']+r"$|$"

    parameter_range = OrderedDict()
    for pn in config.parameter_names:
        parameter_range[pn] = OrderedDict()
        for i in range(n_iterations):
            parameter_range[pn][i] = OrderedDict()
            parameter_range[pn][i]['min'] = kde_data[i].df[pn].min()
            parameter_range[pn][i]['max'] = kde_data[i].df[pn].max()
        parameter_range[pn]['all'] = OrderedDict()
        parameter_range[pn]['all']['min'] = min([parameter_range[pn][i]['min'] for i in range(n_iterations)])
        parameter_range[pn]['all']['max'] = max([parameter_range[pn][i]['max'] for i in range(n_iterations)])
    
    abserr_range = OrderedDict()
    for q in config.qoi_names:
        aen="{}.abserr".format(q)
        abserr_range[aen] = OrderedDict()
        for i in range(n_iterations):
            abserr_range[aen][i] = OrderedDict()
            abserr_range[aen][i]['min'] = kde_data[i].df[aen].min()
            abserr_range[aen][i]['max'] = kde_data[i].df[aen].max()
        abserr_range[aen]['all'] = OrderedDict()
        abserr_range[aen]['all']['min'] = min([abserr_range[aen][i]['min'] for i in range(n_iterations)])
        abserr_range[aen]['all']['max'] = max([abserr_range[aen][i]['max'] for i in range(n_iterations)])

    from pypospack.pyposmat.visualization import Pyposmat2DDensityPlots
    param_2d_kde_plots_dst_dir = os.path.join(
            dst_dir,'param_2d_kde_plots')
    if not os.path.exists(param_2d_kde_plots_dst_dir):
        os.mkdir(param_2d_kde_plots_dst_dir)
        for i in range(n_iterations):
            for j,jpn in enumerate(config.parameter_names):
                for k,kpn in enumerate(config.parameter_names):
                    if j<k:
                        print("kde_2d_plots:{}:{} vs {}".format(i,jpn,kpn))
                        plot_fn = os.path.join(param_2d_kde_plots_dst_dir,"{}__{}__{}.eps".format(jpn,kpn,i))
                        x_lims=[parameter_range[jpn]['all']['min'],parameter_range[jpn]['all']['max']]
                        y_lims=[parameter_range[kpn]['all']['min'],parameter_range[kpn]['all']['max']]
                        plotter=Pyposmat2DDensityPlots()
                        plotter.configuration=config
                        plotter.datafile=kde_data[i]
                        plotter.plot(
                                x_name=jpn,
                                y_name=kpn,
                                x_lims=x_lims,
                                y_lims=y_lims,
                                fn_plot_out=plot_fn)

    abserr_2d_kde_plots_dst_dir = os.path.join(
            dst_dir,'pareto_2d_kde_plots')
    if not os.path.exists(abserr_2d_kde_plots_dst_dir):
        os.mkdir(abserr_2d_kde_plots_dst_dir)
        for i in range(n_iterations):
            for j,jpn in enumerate(config.qoi_names):
                for k,kpn in enumerate(config.qoi_names):
                    if j<k:
                        aen1="{}.abserr".format(jpn)
                        aen2="{}.abserr".format(kpn)
                        print("pareto_2d_plots:{}:{} vs {}".format(i,aen1,aen2))
                        plot_fn = os.path.join(abserr_2d_kde_plots_dst_dir,"{}__{}__{}.eps".format(aen1,aen2,i))
                        x_lims=[abserr_range[aen1]['all']['min'],abserr_range[aen1]['all']['max']]
                        y_lims=[abserr_range[aen2]['all']['min'],abserr_range[aen2]['all']['max']]
                        plotter=Pyposmat2DDensityPlots()
                        plotter.configuration=config
                        plotter.datafile=kde_data[i]
                        plotter.plot(
                                x_name=aen1,
                                y_name=aen2,
                                x_lims=x_lims,
                                y_lims=y_lims,
                                fn_plot_out=plot_fn)
    
                        parameter_range = OrderedDict()
    for pn in config.parameter_names:
        parameter_range[pn] = OrderedDict()
        for i in range(n_iterations):
            parameter_range[pn][i] = OrderedDict()
            parameter_range[pn][i]['min'] = kde_data[i].df[pn].min()
            parameter_range[pn][i]['max'] = kde_data[i].df[pn].max()
        parameter_range[pn]['all'] = OrderedDict()
        parameter_range[pn]['all']['min'] = min([parameter_range[pn][i]['min'] for i in range(n_iterations)])
        parameter_range[pn]['all']['max'] = max([parameter_range[pn][i]['max'] for i in range(n_iterations)])
    
    for pn,pv in parameter_range.items():
        for i,v in pv.items():
            print("{:20} {:5} {:10} {:10}".format(pn,i,v['min'],v['max']))
    
    for qn in config.qoi_names:
        print(qn)

    for qn in config.qoi_names:
        print(qn)

    for i,qn1 in enumerate(config.qoi_names): 
        for j,qn2 in enumerate(config.qoi_names):
            if i<j: print(qn1,qn2)

    for i,pn1 in enumerate(config.parameter_names):
        for j,pn2 in enumerate(config.parameter_names):
            if i<j: print(pn1,pn2)

