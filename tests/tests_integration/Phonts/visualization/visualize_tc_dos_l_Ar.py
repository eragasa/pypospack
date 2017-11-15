import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_spectral_k(filename="tc_dos_l.dat"):
    """
    Reads the spectrial thermal conductivity information
 
    """
    # column headers for the data  
    #tcdosl_labels = [
    #     "wavelength",
    #     "k_xx_raw","k_xx_smooth",
    #     "k_yy_raw","k_yy_smooth",
    #     "k_zz_raw","k_zz_smooth"]

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

def make_tcdosl_plot(
        data_filename='tc_dos_l.dat',
        figure_prefix='tc_dos_l',
        xlim=None,
        ylim=None):
        
    tcdosl_df_dict = read_spectral_k(filename=data_filename)
    for keys in tcdosl_df_dict.keys():
        tcdosl_figure_filename = tcdosl_figure_prefix + '_' + str(keys) + 'K' + '.png'
        figure = plt.figure()
        tcdosl_plot = figure.add_subplot(111)
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_xx_raw'], label='k_xx_raw', color='g')
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_yy_raw'], label='k_yy_raw', color='b')
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_zz_raw'], label='k_zz_raw', color='c')
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_xx_smooth'], label='k_xx_smooth', color='y')
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_yy_smooth'], label='k_yy_smooth', color='m')
        tcdosl_plot.plot(tcdosl_df_dict[keys]['wavelength'],tcdosl_df_dict[keys]['k_zz_smooth'], label='k_zz_smooth', color='r')
        tcdosl_title=plt.title('Spectral thermal conductivity'+ ' at ' + str(keys)+ 'K')
        tcdosl_xlabel=plt.xlabel('Wavelength (Angstrom)')
        tcdosl_ylabel=plt.ylabel('Thermal conductivity (W/mK)')
        tcdosl_legend=plt.legend(loc='upper right', prop={'size':7})
        
        #set axis here
        if xlim is not None:
            tcdosl_plot.set_xlim(xlim)
        if ylim is not None:
            tcdosl_plot.set_ylim(ylim)
        
        figure.savefig(tcdosl_figure_filename)
        plt.close

if __name__ == "__main__":
    phonts_sim_dir = 'Ar_result'
    tcdosl_data_filename = os.path.join(phonts_sim_dir,'tc_dos_l.dat')
    tcdosl_figure_prefix = 'tc_dos_l'    

    assert type(tcdosl_data_filename)
    assert os.path.isfile(tcdosl_data_filename)

    tcdosl_df_dict = read_spectral_k(filename=tcdosl_data_filename)

    #for k,v in tcdosl_df_dict.items():
    #    print(k)
    #    assert type(k) is int
    #    assert type(v) is pd.DataFrame
    #print(tcdosl_df_dict[400]['k_xx_raw'])
    
    # example how to use make_tcdosl_plot(i)
    make_tcdosl_plot(
        data_filename = tcdosl_data_filename,
        figure_prefix = tcdosl_figure_prefix,
        xlim = [0,0.25],
        ylim = [0,2.5])
    


