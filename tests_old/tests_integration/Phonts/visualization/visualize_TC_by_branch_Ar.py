import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_spectral_k(filename="TC_by_branch.dat"):
    """
    Reads the spectrial thermal conductivity information
 
    """
    # column headers for the data  
    #TCbybranch_labels = [
    #     "wavelength",
    #     "k_xx_raw","k_xx_smooth",
    #     "k_yy_raw","k_yy_smooth",
    #     "k_zz_raw","k_zz_smooth",
    #     "lifetime_dos1 ","lifetime_dos2"]

    TCbybranch_labels = [
         "index",
         "k_xx","k_yy","k_zz"]

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
    TCbybranch_dict = OrderedDict()

    for il,line in enumerate(lines):
        if line.startswith('# Temp:'):
            args = line.split(':')
            T = int(float(args[1].strip()))
            temperatures.append(T)
            TCbybranch_dict[T] = subselect_table_block(il,lines)

    TCbybranch_df_dict = OrderedDict()
    for temp in temperatures:
        TCbybranch_df_dict[temp] = pd.DataFrame(
            copy.deepcopy(TCbybranch_dict[temp]),
            columns=list(TCbybranch_labels))

    return {k:v.copy() for k,v in TCbybranch_df_dict.items()}

def make_TCbybranch_plot(
        data_filename='TC_by_branch.dat',
        figure_prefix='TC_by_branch',
        xlim=None,
        ylim=None):
        
    TCbybranch_df_dict = read_spectral_k(filename=data_filename)
    for keys in TCbybranch_df_dict.keys():
        TCbybranch_figure_filename = TCbybranch_figure_prefix + '_' + str(keys) + 'K' + '.png'
        figure = plt.figure()
        TCbybranch_plot = figure.add_subplot(111)
        TCbybranch_plot.plot(TCbybranch_df_dict[keys]['index'],TCbybranch_df_dict[keys]['k_xx'], 'g-o', label='k_xx')
        TCbybranch_plot.plot(TCbybranch_df_dict[keys]['index'],TCbybranch_df_dict[keys]['k_yy'], 'b-^', label='k_yy_raw')
        TCbybranch_plot.plot(TCbybranch_df_dict[keys]['index'],TCbybranch_df_dict[keys]['k_zz'], 'c-s', label='k_zz_raw')
        TCbybranch_title=plt.title('Spectral thermal conductivity'+ ' at ' + str(keys)+ 'K')
        TCbybranch_xlabel=plt.xlabel('Branch Number')
        TCbybranch_ylabel=plt.ylabel('Thermal conductivity (W/mK)')
        TCbybranch_legend=plt.legend(loc='upper right', prop={'size':7})
        
        #set axis here
        if xlim is not None:
            TCbybranch_plot.set_xlim(xlim)
        if ylim is not None:
            TCbybranch_plot.set_ylim(ylim)
        
        figure.savefig(TCbybranch_figure_filename)
        plt.close

if __name__ == "__main__":
    phonts_sim_dir = 'Ar_result'
    TCbybranch_data_filename = os.path.join(phonts_sim_dir,'TC_by_branch.dat')
    TCbybranch_figure_prefix = 'TC_by_branch'    

    assert type(TCbybranch_data_filename)
    assert os.path.isfile(TCbybranch_data_filename)

    TCbybranch_df_dict = read_spectral_k(filename=TCbybranch_data_filename)

    #for k,v in TCbybranch_df_dict.items():
    #    print(k)
    #    assert type(k) is int
    #    assert type(v) is pd.DataFrame
    #print(TCbybranch_df_dict[400]['k_xx_raw'])
    
    # example how to use make_TCbybranch_plot(i)
    make_TCbybranch_plot(
        data_filename = TCbybranch_data_filename,
        figure_prefix = TCbybranch_figure_prefix,
        xlim = [0,14],
        ylim = [0,7])
    


