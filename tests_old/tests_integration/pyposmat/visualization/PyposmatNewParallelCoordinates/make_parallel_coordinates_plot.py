import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile, PyposmatConfigurationFile
from pypospack.pyposmat.visualization.parallel_plot_new import PyposmatParallelCoordinatesPlot


if __name__ == "__main__":
    ppcp = PyposmatParallelCoordinatesPlot()
    configfile = PyposmatConfigurationFile()
    # sorry about the absolute paths
    configfile.read("/home/seaton/python-repos/pypospack/data/Si__sw__data/pareto_optimization_p_3.5_q_0.5/pyposmat.config.in")
    datafile = PyposmatDataFile()
    datafile.read("/home/seaton/python-repos/pypospack/data/Si__sw__data/pareto_optimization_p_3.5_q_0.5/pyposmat.kde.5.out")
    datafile.create_normalized_errors(qoi_targets=configfile.qoi_targets)
    df = datafile.df
    ppcp.add_dataframe("blue", "kde5", df, names=datafile.normalized_names)
    datafile.read("/home/seaton/python-repos/pypospack/data/Si__sw__data/pareto_optimization_p_3.5_q_0.5/pyposmat.kde.15.out")    
    datafile.create_normalized_errors(qoi_targets=configfile.qoi_targets)
    ppcp.add_datafile("orange", "kde15", datafile, names=datafile.normalized_names)
    ppcp.make_plot(filename="parallel_plot.png", xlabels=datafile.normalized_names, 
                   ylabel="% error", title="Si sw", ylim=(-175, 25), legend_loc="lower right")
