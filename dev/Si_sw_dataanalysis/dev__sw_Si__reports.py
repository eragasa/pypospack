import os,shutil,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatDataAnalyzer

def make_histogram(data,dst_file):
    plt.close('all')
    f,ax = plt.subplots()
    ax.hist(data,normed=1)
    f.savefig(dst_file)

if __name__ == "__main__":
    import Si_sw

#WRITE CONFIGURATION FILE

    Si_sw_configuration = PyposmatConfigurationFile()
    Si_sw_configuration.qois = Si_sw.Si_sw_qoi_db.qois
    Si_sw_configuration.potential = Si_sw.potential
    Si_sw_configuration.structures = Si_sw.Si_sw_structures
    Si_sw_configuration.sampling_type = Si_sw.Si_sw_sampling
    Si_sw_configuration.sampling_distribution = Si_sw.Si_sw_param_dist
    Si_sw_configuration.write(filename = 'pypospack.config.in')
    Si_sw_configuration.read(filename = 'pypospack.config.in')

    pypospack_filename_in = 'pypospack.config.in'

    datafile_filename = os.path.join('data','pyposmat.kde.2.out')
    datafile = PyposmatDataFile(filename=datafile_filename)
    datafile.read()

    histogram_dir = 'histograms'
    if os.path.isdir(histogram_dir):
        shutil.rmtree(histogram_dir)
    os.mkdir(histogram_dir)

    for p_name in datafile.parameter_names:
        _data = datafile.df[p_name],
        _filename = os.path.join(
                histogram_dir,
                '{}.png'.format(p_name)
            )
       
        make_histogram(
                data=_data,
                dst_file=_filename)

    for pidx1,pn1 in enumerate(datafile.parameter_names):
        for pidx2,pn2 in enumerate(datafile.parameter_names):
            if pidx1 < pidx2:
                _data1 = datafile.df[pn1]
                _data2 = datafile.df[pn2]
    #qoi_names
    #error_names
