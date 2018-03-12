# imports here
'''
TODO:
scaling and box zoom

probability density plot
-histogram
1d fit
2s fit
-multivariant normal distribution

kernel density estimate
-pyflames post module

principal components analysis
clustering
manifold learning T-sne
'''
import os
from pypospack.pyposmat.visualization.bokeh import PyposmatBokehVisualizer

if __name__ == "__main__":
    filename = os.path.join(
            'data__Ni__eam__born_exp_fs_02',
            'pyposmat.kde.9.out')
    #filename = os.path.join(
    #        'data__Ni__eam__born_exp_fs_03',
    #        "pyposmat.kde.0.out")

    vizdemo = PyposmatBokehVisualizer()
    vizdemo.load_data_file(filename= filename)
    print(vizdemo.err_names[0])
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
