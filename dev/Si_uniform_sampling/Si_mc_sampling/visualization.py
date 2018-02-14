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
import pypospack.visualization as visualization

if __name__ == "__main__":
    data_dir = 'Si_uniform_sampling'
    filename = 'pypospack.results.out'

    print(type(data_dir))
    print(type(filename))
    print(type(os.path.join(data_dir,filename)))
    vizdemo = visualization.ParetoOptimizationParamVsErrorScatter()
    vizdemo.load_data_file(fname= filename)
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
