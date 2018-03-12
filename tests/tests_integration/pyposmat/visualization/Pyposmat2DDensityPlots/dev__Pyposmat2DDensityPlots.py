import os
from collections import OrderedDict

from  pypospack.pyposmat.visualization import Pyposmat2DDensityPlots

if __name__ == "__main__":
    data_dir = os.path.join(
            "../../../../",
            "data_test",
            "Ni__eam__born_exp_fs_00",
            "data__Ni__eam__born_exp_fs_02"
        )
    fn_config = os.path.join(data_dir,"pyposmat.config.in")
    fn_results = os.path.join(data_dir,"pyposmat.kde.9.out")
    myplot = Pyposmat2DDensityPlots()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.plot(
        x_name='p_NiNi_r0',
        y_name='d_Ni_r0',
        XY_density_h='chiu1999',
    )
