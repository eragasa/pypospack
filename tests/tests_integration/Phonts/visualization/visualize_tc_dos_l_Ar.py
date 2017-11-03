import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import phonts
from collections import OrderedDict

def read_spectral_k(filename="tc_dos_l.dat"):
    """
    Reads the spectrial thermal conductivity information
 
    """
    # column headers for the data  
    tcdosl_labels = [
         "wavelength",
         "k_xx_raw","k_yy_raw","k_zz_raw",
         "k_xx_smooth","k_yy_smooth","k_zz_smooth"]

    def subselect_table_block(i_start,lines):
        i = i_start + 1

        table = []
        while(lines[i].strip() != ""):
            args = lines[i].split()
            args = [arg.strip() for arg in args]
            args = [float(arg) for arg in args]
            table.append(args)
            i += 1  
        return np.array(table)

    line = None # initialize
    with open(filename,'r') as f:
        lines = f.readlines()
    lines = [s.strip() for s in lines]

    temperatures = []
    tcdosl_dict = OrderedDict()

    for il,line in enumerate(lines):
        if line.startswith('# Temp:'):
            args = line.split(':')
            T = int(float(args[1].strip()))
            temperatures.append(T)
            tcdosl_dict[T] = subselect_table_block(il,lines)

    tcdosl_df_dict = OrderedDict()
    for temp in temperatures:
        tcdosl_df_dict[temp] = pd.DataFrame(
            copy.deepcopy(tcdosl_dict[temp]),
            columns=list(tcdosl_labels))

    return {k:v.copy() for k,v in tcdosl_df_dict.items()}

if __name__ == "__main__":
    phonts_sim_dir = 'Ar_result'
    tcdosl_filename = os.path.join(phonts_sim_dir,'tc_dos_l.dat')

    assert type(tcdosl_filename)
    assert os.path.isfile(tcdosl_filename)

    tcdosl_df_dict = read_spectral_k(filename=tcdosl_filename)

    for k,v in tcdosl_df_dict.items():
        print(k)
        assert type(k) is int
        assert type(v) is pd.DataFrame
    print(tcdosl_df_dict[400]['wavelength'])
