import os,shutil
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.visualization import Pyposmat2DDensityPlot

data_directory = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','Si__sw__data',
        'pareto_optimization_unconstrained'
        )

plot_directory = 'qoi_2d_density_plots'
if os.path.isdir(plot_directory):
    shutil.rmtree(plot_directory)
os.mkdir(plot_directory)

o_config = PyposmatConfigurationFile()
o_config.read(filename=os.path.join(data_directory,'pyposmat.config.in'))

max_iteration = o_config.n_iterations
o_data = PyposmatDataFile()
o_data.read(filename=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(max_iteration)))

for i, name_1 in enumerate(o_config.qoi_names):
    for j, name_2 in enumerate(o_config.qoi_names):
        if i < j:
            plot_fn = '{}__{}'.format(name_1,name_2).replace('.','_')
            plot_fn += ".png"
            plot_fn = os.path.join(plot_directory,plot_fn)
            o = Pyposmat2DDensityPlot(config=o_config,data=o_data)
            o.plot(x_name=name_1,y_name=name_2)
            o.fig.savefig(plot_fn,dpi=1300)
