import os
from pypospack.pyposmat.visualization.plots_1d import Pyposmat1DHistogramWithDensityPlots

pyposmat_root = [v for v in os.environ['PYTHONPATH'].strip().split(':') if v.endswith('pypospack')][0]
data_directory = os.path.join(pyposmat_root,'data','MgO_pareto_data')

config_fn = os.path.join(data_directory,'pyposmat.config.in')
data_fn = os.path.join(data_directory,"culled_004.out")


myplot = Pyposmat1DHistogramWithDensityPlots()
myplot.read_configuration(filename=config_fn)
myplot.read_datafile(filename=data_fn)

plot_fn = "sum_abs_error.eps"
myplot.plot(
        x_name="sum_all.nerr",
        include_histogram=True,
        include_kde=True,
        include_normal=False,
        filename=plot_fn)

plot_fn = "sum_sq_error.eps"
myplot.plot(
        x_name="sum_sq.nerr",
        include_histogram=True,
        include_kde=True,
        include_normal=False,
        filename=plot_fn)
