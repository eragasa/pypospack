import pytest
import os
from pypospack.pyposmat.visualization.rugplot_pareto import PyposmatParetoRugplot

pypospack_root_dir = [v for v in os.environ['PYTHONPATH'].split(':') if v.endswith('pypospack')][0]
# -----------------------------------------------------------------------------
# DEFINE WHERE TO FIND ALL THE DATA
# -----------------------------------------------------------------------------
datafile_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/qoiplus_005.out')
config_fn = os.path.join(pypospack_root_dir,'examples/MgO__buck__add_additional_qoi/data/pyposmat.config.in')

# -----------------------------------------------------------------------------
# DEFINE WHERE TO PUT ALL THE OUTPUT
# -----------------------------------------------------------------------------
output_directory = "./"
output_plot_fn = os.path.join(output_directory,'rugplot_MgO_buck.png')

if __name__ == "__main__":

    o_rugplot = PyposmatParetoRugplot()
    o_rugplot.read_configuration(filename=config_fn)
    o_rugplot.read_datafile(filename=datafile_fn)
    o_rugplot.make_plot(
            filename=output_plot_fn,
            include_qois=True,
            include_qois_v=True,
            qoi_excluded_names=None)
    exit()
