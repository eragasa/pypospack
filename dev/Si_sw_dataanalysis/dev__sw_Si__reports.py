import os,shutil,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatDataAnalyzer

if __name__ == "__main__":
    import Si_sw

#WRITE CONFIGURATION FILE

    #Si_sw_configuration = PyposmatConfigurationFile()
    #Si_sw_configuration.qois = Si_sw.Si_sw_qoi_db.qois
    #Si_sw_configuration.potential = Si_sw.potential
    #Si_sw_configuration.structures = Si_sw.Si_sw_structures
    #Si_sw_configuration.sampling_type = Si_sw.Si_sw_sampling
    #Si_sw_configuration.sampling_distribution = Si_sw.Si_sw_param_dist
    #Si_sw_configuration.write(filename = 'pypospack.config.in')
    #Si_sw_configuration.read(filename = 'pypospack.config.in')

    pypospack_filename_in = 'pypospack.config.in'

    #datafile_filename = os.path.join('data','pyposmat.kde.2.out')
    #datafile = PyposmatDataFile(filename=datafile_filename)
    #datafile.read()

    #print(type(datafile.df))
    #print(datafile.df)

    plt.close('all')
    f,ax = plt.subplots()

    #p_name = 'SiSiSi_epsilon'
    #ax[0].hist(datafile.df[p_name],normed=1)
    
    #for p_name in Si_sw_configuration.parameter_names:
    #    print(p_name)
