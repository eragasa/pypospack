import os
import numpy as np
import matplotlib.pyplot as plt

def get_data_from_phonts_file(filename):
    def process_first_line(line):
        args = line.strip().split()
        args = [arg.strip() for arg in args]
        args = args[1:]
        return args

    labels = None
    data = None
    with open(filename) as f:
        lines = f.readlines()
    # except FileNotFoundError

   #initialize variables
    values = []
    for i,line in enumerate(lines):
        if i == 0:
            labels = process_first_line(line)
        else:
            args = line.strip().split()
            values.append([float(arg) for arg in args])

    data = np.array(values)

    return data,labels

def get_freq_dispersion_data(filename='freq_dispersion.dat'):
     (freq_dispersion_data,freq_dispersion_labels) = get_data_from_phonts_file(filename)
     return freq_dispersion_data,freq_dispersion_labels

def make_freq_dispersion_plot(
        data_filename='freq_dispersion.dat',
        figure_filename='freq_dispersion.png',
        xlim=None,
        ylim=None):
    """
    a short descript here 

    a long description here. alskdjflakdjflkjldkajflakjflkadjflkjlfkdjlkjkjlkj
    aldskfjlkdjflkjljkalfdjlkdjflaskjdfljkdfl.

    Args:
        data_filename(str): input file from PhonTS which has density of state
            (DOS) data on it.  By default is it is set to 'freq_dispersion.dat'
        figure_filename(str): output file to save plot to
        xlim(list of float): Controls the range of the x-axis. [xmin,xmax]. 
            By default this is set to None, and the range is determined 
            automatically
        ylim(list of float): Controls the range of the y-axis. [ymin,ymax].
            By default this is set to None, and the rnage is determined
            automatically.

    Returns:
        None

    """
    (freq_dispersion_data,freq_dispersion_labels) = get_freq_dispersion_data(filename=data_filename)

    figure = plt.figure()
    freq_dispersion_plot = figure.add_subplot(111)
    
    for i in range(4,16):
        freq_dispersion_plot.plot(freq_dispersion_data[:,0],freq_dispersion_data[:,i])
    
    freq_dispersion_title=plt.title('Calculated phonon dispersion relation of Argon', fontname='Times New Roman')
    freq_dispersion_ylabel=plt.ylabel('Frequency (meV)', fontname='Times New Roman')
    freq_dispersion_xticks2=plt.xticks([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 270, 300],[r'Gamma', r'X', r'W', r'K', r'Gamma', r'L', r'U', r'W', r'L', r'K', r'U', r'X'], fontname='Times New Roman')
    
    for i in range(1,11):
        freq_dispersion_vlines=plt.vlines(i*30, 0, 12)    
    
    freq_dispersion_font=plt.rc('font', family='Times New Roman')
    # set axis here
    if xlim is not None:
        freq_dispersion_plot.set_xlim(xlim)
    if ylim is not None:
        freq_dispersion_plot.set_ylim(ylim)
     
    figure.savefig(figure_filename)
    plt.close()
    
if __name__ == "__main__":

    phonts_sim_dir = 'Ar_result'
    freq_dispersion_data_filename = os.path.join(phonts_sim_dir,'freq_dispersion.dat')
    freq_dispersion_figure_filename = 'freq_dispersion.png'

    assert type(freq_dispersion_data_filename) is str
    assert os.path.isfile(freq_dispersion_data_filename)

    # example how to use get_freq_dispersion_data()
    #(freq_dispersion_data,freq_dispersion_labels) = get_freq_dispersion_data(
    #    filename=freq_dispersion_data_filename)

    # example how to use make_freq_dispersion_plot(i)
    make_freq_dispersion_plot(
        data_filename = freq_dispersion_data_filename,
        figure_filename = freq_dispersion_figure_filename,
        xlim=[0,300],
        ylim=[0,12])


