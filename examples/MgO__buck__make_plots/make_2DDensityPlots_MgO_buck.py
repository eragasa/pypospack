import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace

import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from matplotlib import rc

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization
from pypospack.pyposmat.visualization import Pyposmat2DDensityPlots
from scipy.stats import gaussian_kde

data_directory = "../../data/MgO_pareto_data"
fn_config = os.path.join(data_directory,'pyposmat.config.in')
fn_results = os.path.join(data_directory,"culled_005.out")

myplot = Pyposmat2DDensityPlots()
myplot.read_datafile(filename=fn_results)
myplot.read_configuration(filename=fn_config)
myplot.plot(
    x_name='chrg_Mg',
    y_name='MgO_A'
)
