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
  
    error_names = myplot.configuration.error_names
    abserror_names = ['{}.abserr'.format(n) for n in error_names]
    
    for ie,ien in enumerate(abserror_names):
        for je,jen in enumerate(abserror_names):
            if ie < je:
                x_name=ien
                y_name_jen

                print(ie,je,ien,jen)
    #myplot.plot(
    #
    #    x_name='Ni_fcc.G.abserr',
    #    y_name='Ni_fcc.B.abserr'
    #)
