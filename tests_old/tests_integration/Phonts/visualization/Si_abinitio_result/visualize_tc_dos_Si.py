import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_spectral_k(filename="tc_dos.dat"):
    """
    Reads the spectrial thermal conductivity information
 
    """

    tcdos_labels = [
         "wavelength",
         "k_xx_raw","k_yy_raw","k_zz_raw",
         "k_xx_smooth","k_yy_smooth","k_zz_smooth",
         "lifetime_dos1 ","lifetime_dos2"]

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
    tcdos_dict = OrderedDict()

    for il,line in enumerate(lines):
        if line.startswith('# Temp:'):
            args = line.split(':')
            T = int(float(args[1].strip()))
            temperatures.append(T)
            tcdos_dict[T] = subselect_table_block(il,lines)

    tcdos_df_dict = OrderedDict()
    for temp in temperatures:
        tcdos_df_dict[temp] = pd.DataFrame(
            copy.deepcopy(tcdos_dict[temp]),
            columns=list(tcdos_labels))

    return {k:v.copy() for k,v in tcdos_df_dict.items()}

def normalize_tcdos(
        data_filename='tc_dos.dat'):

    tcdos_df_dict = read_spectral_k(filename=data_filename)
    tcdos_df_dict_n = tcdos_df_dict

    for k, v in tcdos_df_dict.items():

        k_xx_raw = sum(list(tcdos_df_dict[k]['k_xx_raw']))
        k_yy_raw = sum(list(tcdos_df_dict[k]['k_yy_raw']))
        k_zz_raw = sum(list(tcdos_df_dict[k]['k_zz_raw']))
        k_xx_smooth = sum(list(tcdos_df_dict[k]['k_xx_smooth']))
        k_yy_smooth = sum(list(tcdos_df_dict[k]['k_yy_smooth']))
        k_zz_smooth = sum(list(tcdos_df_dict[k]['k_zz_smooth']))
        
        tcdos_df_dict_n[k]['k_xx_raw'] = tcdos_df_dict[k]['k_xx_raw']/k_xx_raw
        tcdos_df_dict_n[k]['k_yy_raw'] = tcdos_df_dict[k]['k_yy_raw']/k_yy_raw
        tcdos_df_dict_n[k]['k_zz_raw'] = tcdos_df_dict[k]['k_zz_raw']/k_zz_raw
        tcdos_df_dict_n[k]['k_xx_smooth'] = tcdos_df_dict[k]['k_xx_smooth']/k_xx_smooth
        tcdos_df_dict_n[k]['k_yy_smooth'] = tcdos_df_dict[k]['k_yy_smooth']/k_yy_smooth
        tcdos_df_dict_n[k]['k_zz_smooth'] = tcdos_df_dict[k]['k_zz_smooth']/k_zz_smooth

    return {k:v.copy() for k,v in tcdos_df_dict_n.items()}

def make_tcdos_plot(
        data_filename='tc_dos.dat',
        figure_prefix='tc_dos',
        xlim=None,
        ylim=None):
        
    tcdos_df_dict = normalize_tcdos(data_filename=data_filename)
    for keys in tcdos_df_dict.keys():
        tcdos_figure_filename = tcdos_figure_prefix + '_' + str(keys) + 'K' + '.png'
        figure = plt.figure()
        tcdos_plot = figure.add_subplot(111)
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_xx_raw'], label='k_xx_raw', color='g')
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_yy_raw'], label='k_yy_raw', color='b')
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_zz_raw'], label='k_zz_raw', color='c')
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_xx_smooth'], label='k_xx_smooth', color='y')
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_yy_smooth'], label='k_yy_smooth', color='m')
        tcdos_plot.plot(tcdos_df_dict[keys]['wavelength'],tcdos_df_dict[keys]['k_zz_smooth'], label='k_zz_smooth', color='r')
        tcdos_title=plt.title('Spectral thermal conductivity'+ ' at ' + str(keys)+ 'K', fontname='Times New Roman')
        tcdos_xlabel=plt.xlabel('Frequency (THz)', fontname='Times New Roman')
        tcdos_ylabel=plt.ylabel('Thermal conductivity (W/mK)', fontname='Times New Roman')
        tcdos_legend=plt.legend(loc='upper right', prop={'size':8})
        
        tcdos_font=plt.rc('font', family='Times New Roman')
        #set axis here
        if xlim is not None:
            tcdos_plot.set_xlim(xlim)
        if ylim is not None:
            tcdos_plot.set_ylim(ylim)
        
        figure.savefig(tcdos_figure_filename)
        plt.close

if __name__ == "__main__":
    phonts_sim_dir = 'Si_abinitio_result'
    tcdos_data_filename = os.path.join(phonts_sim_dir,'tc_dos.dat')
    tcdos_figure_prefix = 'tc_dos'    

    assert type(tcdos_data_filename)
    assert os.path.isfile(tcdos_data_filename)

    tcdos_df_dict = read_spectral_k(filename=tcdos_data_filename)
    
    # example how to use make_tcdos_plot(i)
    make_tcdos_plot(
        data_filename = tcdos_data_filename,
        figure_prefix = tcdos_figure_prefix,
        xlim = [0,18],
        ylim = [0,0.1])
    


