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

def get_pdos_data(filename='pdos.dat'):
     (pdos_data,pdos_labels) = get_data_from_phonts_file(filename)
     return pdos_data,pdos_labels

def make_pdos_plot(
        data_filename='pdos.dat',
        figure_filename='pdos.png',
        xlim=None,
        ylim=None):
    """
    a short descript here 

    a long description here. alskdjflakdjflkjldkajflakjflkadjflkjlfkdjlkjkjlkj
    aldskfjlkdjflkjljkalfdjlkdjflaskjdfljkdfl.

    Args:
        data_filename(str): input file from PhonTS which has density of state
            (DOS) data on it.  By default is it is set to 'pdos.dat'
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
    (pdos_data,pdos_labels) = get_pdos_data(filename=data_filename)

    figure = plt.figure()
    pdos_plot = figure.add_subplot(111)
    pdos_plot.plot(pdos_data[:,2],pdos_data[:,3], label='density of states as calculated')
    pdos_plot.plot(pdos_data[:,2],pdos_data[:,4], label='density of states as smoothed')
    pdos_title=plt.title('Density of states for Argon')
    pdos_xlabel=plt.xlabel('Frequency (meV)')
    pdos_ylabel=plt.ylabel('PDOS')
    pdos_legend=plt.legend(loc='upper right', prop={'size':7})

    # set axis here
    if xlim is not None:
        pdos_plot.set_xlim(xlim)
    if ylim is not None:
        pdos_plot.set_ylim(ylim)
     
    figure.savefig(figure_filename)
    plt.close()

if __name__ == "__main__":

    phonts_sim_dir = 'Ar_result'
    pdos_data_filename = os.path.join(phonts_sim_dir,'pdos.dat')
    pdos_figure_filename = 'pdos.png'

    assert type(pdos_data_filename) is str
    assert os.path.isfile(pdos_data_filename)

    # example how to use get_pdos_data()
    (pdos_data,pdos_labels) = get_pdos_data(
        filename=pdos_data_filename)

    # example how to use make_pdos_plot(i)
    make_pdos_plot(
        data_filename = pdos_data_filename,
        figure_filename = pdos_figure_filename,
        xlim=[0,600],
        ylim=None)


