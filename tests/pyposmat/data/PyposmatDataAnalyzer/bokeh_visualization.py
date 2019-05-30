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
    data_dir = 'data'
    filename = os.path.join("resources",'pyposmat.results.0a.out')

    vizdemo = PyposmatBokehVisualizer()
    vizdemo.load_data_file(filename= filename)
    print(vizdemo.err_names[0])
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
