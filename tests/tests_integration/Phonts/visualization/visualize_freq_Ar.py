import os
import numpy as np
import matplotlib.pyplot as plt

# Haven't been finished
# need to use the data to obtain the figure of dispersion relation

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

def get_freq_data(filename='freq.dat'):
    (freq_data,freq_labels) = get_data_from_phonts_file(filename)
    return freq_data,freq_labels

def make_freq_plot(
        data_filename='freq.dat',
        figure_filename='freq.png',
        xlim=None,
        ylim=None):
    """
    a short descript here 

    a long description here. alskdjflakdjflkjldkajflakjflkadjflkjlfkdjlkjkjlkj
    aldskfjlkdjflkjljkalfdjlkdjflaskjdfljkdfl.

    Args:
        data_filename(str): input file from PhonTS which has density of state
            (DOS) data on it.  By default is it is set to 'freq.dat'
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
    (freq_data,freq_labels) = get_freq_data(filename=data_filename)

    figure = plt.figure()
    freq_plot = figure.add_subplot(111)
    freq_plot.plot(freq_data[:,0],freq_data[:,1], label='frequency 1')
    freq_plot.plot(freq_data[:,0],freq_data[:,2], label='frequency 2')
    freq_plot.plot(freq_data[:,0],freq_data[:,3], label='frequency 3')
    freq_plot.plot(freq_data[:,0],freq_data[:,4], label='frequency 4')
    freq_plot.plot(freq_data[:,0],freq_data[:,5], label='frequency 5')
    freq_plot.plot(freq_data[:,0],freq_data[:,6], label='frequency 6')
    freq_plot.plot(freq_data[:,0],freq_data[:,7], label='frequency 7')
    freq_plot.plot(freq_data[:,0],freq_data[:,8], label='frequency 8')
    freq_plot.plot(freq_data[:,0],freq_data[:,9], label='frequency 9')
    freq_plot.plot(freq_data[:,0],freq_data[:,10], label='frequency 10')
    freq_plot.plot(freq_data[:,0],freq_data[:,11], label='frequency 11')
    freq_plot.plot(freq_data[:,0],freq_data[:,12], label='frequency 12')

    plt.axis(auto)
    #freq_title=plt.title('Dispersion relation for Argon')
    #freq_xlabel=plt.xlabel('Wavevector')
    #freq_ylabel=plt.ylabel('Frequency (THz)')
    freq_legend=plt.legend(loc='upper right', prop={'size':7})

    # set axis here
    if xlim is not None:
        freq_plot.set_xlim(xlim)
    if ylim is not None:
        freq_plot.set_ylim(ylim)
     
    figure.savefig(figure_filename)
    plt.close()

if __name__ == "__main__":

    phonts_sim_dir = 'Ar_result'
    freq_data_filename = os.path.join(phonts_sim_dir,'freq.dat')
    freq_figure_filename = 'freq.png'

    assert type(freq_data_filename) is str
    assert os.path.isfile(freq_data_filename)

    # example how to use get_freq_data()
    (freq_data,freq_labels) = get_freq_data(
        filename=freq_data_filename)

    # example how to use make_freq_plot(i)
    make_freq_plot(
        data_filename = freq_data_filename,
        figure_filename = freq_figure_filename,
        xlim=[0,600],
        ylim=None)


