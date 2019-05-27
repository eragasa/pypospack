import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import pypospack.io.phonts as phonts

if __name__ == "__main__":
    phonts_sim_dir = 'Ar_result'
    freq_data_filename = os.path.join(
             phonts_sim_dir,
             'freq.dat')
    phon_lifetime_data_filename = os.path.join(
             phonts_sim_dir,
             'phon_lifetime.dat')

    bte_data = phonts.PhontsBteData(natoms=4,directory=phonts_sim_dir)
    bte_data.read()
    bte_data.build_data_at_temp(temp=400)
    ph_freq = bte_data.data[400][:,4]
    ph_lt = bte_data.data[400][:,5]
    
    idx_not_zero = np.where(ph_lt != 0)[0]
    ph_freq = bte_data.data[400][idx_not_zero,4]
    ph_inv_lt = 1/bte_data.data[400][idx_not_zero,5]
