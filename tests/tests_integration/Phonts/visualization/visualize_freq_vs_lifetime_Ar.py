import os,copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

class PhontsBteData(object):
    """
    Args:
        directory(str): the directory where the phonts simulation
            output files exist
        natoms(int):number of atoms in the simulation cel
    Attributes:
        directory(str): the directory where the phonts simulation
            output files exist
        ph_lt_file(pypospack.io.phonts.PhononLifetimeFile)
        ph_freq_file(pypospack.io.phonts.PhononFrequencyFile)
        natoms(int):number of atoms in the simulation cel
    """

    def __init__(self,directory,natoms):
        self.directory =directory
        self.ph_lt_file = None
        self.ph_freq_file = None

        self.ph_lt_filename = 'phon_lifetime.dat'
        self.ph_freq_filename = 'freq.dat'

        self.natoms = natoms
        self.data = None

    def read(self,directory=None):
        if directory is not None:
            self.directory = directory
            
        self.ph_lt_file = PhononLifetimeFile(
            natoms=self.natoms,
            filename=os.path.join(
                self.directory,
                self.ph_lt_filename))
        self.ph_freq_file = PhononFrequencyFile(
            natoms=self.natoms,
            filename=os.path.join(
                self.directory,
                self.ph_freq_filename))
        
        self.ph_lt_file.read()
        self.ph_freq_file.read()

    def build_data(self):
        self.data = OrderedDict()
        for temp in self.ph_lt_file.temp:
            print('temp={}'.format(temp))
            self.data[temp] = None
            self.build_data_at_temp(temp)


    def build_data_at_temp(self,temp):
        """
        This method consolidates the data contained in two different files
        so that we can compare phonon frequencies with with phonon lifetimes.
        
        Args:
            temp(int): the temperature from the BTE calculation, this value
                must be in list of values contained in the list ph_lt_file.temp
        """
        if self.data is None:
            self.data = OrderedDict()

        self.data[temp] = [] # initialize our list
        
        self.kpoint_keys_format = "{kp1:.6f}_{kp2:.6f}_{kp3:.6f}"

        freq_n_rows, freq_n_cols = self.ph_freq_file.data.shape
        lt_n_rows, lt_n_cols = self.ph_lt_file.data[temp].shape

        #these indices are the column index for the columns kp1,kp2,kp3 in 
        #phonon frequency file (ph_fr_file) 
        freq_kp1_idx = self.ph_freq_file.col_names.index('kp1')
        freq_kp2_idx = self.ph_freq_file.col_names.index('kp2')
        freq_kp3_idx = self.ph_freq_file.col_names.index('kp3')

        # these indices are the the column indices for the columns kp1,kp2,kp3
        # in the lifetime frequency file (ph_lt_file)
        lt_kp1_idx = self.ph_lt_file.col_names.index('kp1')
        lt_kp2_idx = self.ph_lt_file.col_names.index('kp2')
        lt_kp3_idx = self.ph_lt_file.col_names.index('kp3')

        for i in range(freq_n_rows):
           freq_kpoint_key = self.kpoint_keys_format.format(
               kp1 = self.ph_freq_file.data[i,freq_kp1_idx],
               kp2 = self.ph_freq_file.data[i,freq_kp2_idx],
               kp3 = self.ph_freq_file.data[i,freq_kp3_idx])
           for j in range(lt_n_rows):
               lt_kpoint_key = self.kpoint_keys_format.format(
                   kp1 = self.ph_lt_file\
                          .data[temp][j,lt_kp1_idx],
                   kp2 = self.ph_lt_file\
                          .data[temp][j,lt_kp2_idx],
                   kp3 = self.ph_lt_file\
                          .data[temp][j,lt_kp3_idx])

               if freq_kpoint_key == lt_kpoint_key:
                   for k in range(3*self.natoms):
                       # here we are building the row for the phonon frequency
                       # and the associated limetime with that phonon
                       # ph_id(int) - unique integer assigned to a phonon for
                       #        identification
                       # kp1,kp2,kp3 - the location of the kpoint associated with
                       #     phonon frequency represented in the basis of the 
                       #     reciprocal lattice
                       # fr - this is the frequency of the phonon in meV
                       # lt - this is the phonon lifetime in ps
                       ph_id = len(self.data[temp])
                       kp1 = self.ph_lt_file\
                           .data[temp][j,lt_kp1_idx]
                       kp2 = self.ph_lt_file\
                           .data[temp][j,lt_kp2_idx]
                       kp3 = self.ph_lt_file\
                           .data[temp][j,lt_kp3_idx]

                       # we need the index associated with the phonon freq
                       fr_idx = self.ph_freq_file.col_names.index(
                           "freq{}".format(k+1))

                       # wew need the index associated with the phonon lifetime
                       lt_idx = self.ph_lt_file.col_names.index(
                           "lt{}".format(k+1))
                       fr = self.ph_freq_file.data[i,fr_idx]
                       lt = self.ph_lt_file.data[temp][j,lt_idx]
                       self.data[temp].append([
                           ph_id,
                           kp1,kp2,kp3,
                           fr,lt])
        self.data[temp] = np.array(self.data[temp])
class PhononLifetimeFile(object):
    
    def __init__(self,natoms,filename='phon_lifetime.dat'):
       self.col_names = ['index','kp1','kp2','kp3']\
                + ['lt{}'.format(i+1) for i in range(3*natoms)]
       self.natoms = natoms
       self.filename = filename
       self.data = read_phon_lifetime(self.filename)
       self.temp = [k for k,v in self.data.items()]

    def read(self):
       self.data = read_phon_lifetime(self.filename)
    
    def print(self):
       for k,v in self.data.items():
           print(k,v.shape)

class PhononFrequencyFile(object):
    
    def __init__(self,natoms,filename='freq.dat'):
        self.col_names = ['index','kp1','kp2','kp3']\
                + ['freq{}'.format(i+1) for i in range(3*natoms)]
        self.filename = filename
        self.data = None
        self.natoms = natoms

    def read(self,filename=None):
        if filename is not None:
            self.filename = filename

        try:
            with open(self.filename,'r') as f:
                lines = f.readlines()
        except FileNotFoundErr as e:
            raise
  
        values_all = []
        for i,line in enumerate(lines):
            args = line.strip().split()
            values_all.append([float(arg) for arg in args])
   
        self.data = np.array(values_all)

    def print(self):
        print(self.data.shape)

def get_data_from_phonts_file(filename):
    #def process_first_line(line):
    #    args = line.strip().split()
    #    args = [arg.strip() for arg in args]
    #    args = args[1:]
    #    return args
    
    #labels = None
    data = None
    data_all = None
    with open(filename) as f:
        lines = f.readlines()
    # except FileNotFoundErr
  
    #initialize variables
    values_all = []
    for i,line in enumerate(lines):
       # if i == 0:
           # labels = process_first_line(line)
       # else:
           args = line.strip().split()
           values_all.append([float(arg) for arg in args])
   
    data_all = np.array(values_all)
    
    return data_all

def get_freq_data(filename='freq.dat'):
    freq_data = get_data_from_phonts_file(filename)
    return freq_data


def read_phon_lifetime(filename='phon_lifetime.dat'):
    """
    Reads phonon lifetime information
    
    """

    def subselect_table_block(i_start,lines):
        i=i_start+1
        
        table = []
        while(lines[i].strip() !=""):
            args = lines[i].split()
            args = [arg.strip() for arg in args]
            args = [float(arg) for arg in args]
            table.append(args)
            i += 1
        return np.array(table)
    
    line = None #initialize
    with open(filename,'r') as f:
        lines = f.readlines()
    lines = [s.strip() for s in lines]
    
    temperatures = []
    phon_lifetime = OrderedDict()

    for il,line in enumerate(lines):
        if line.startswith('# Temp:'):
            args = line.split(':')
            T = int(float(args[1].strip()))
            temperatures.append(T)
            phon_lifetime[T] = subselect_table_block(il,lines)

    return {k:v.copy() for k,v in phon_lifetime.items()}

if __name__ == "__main__":
    phonts_sim_dir = 'Ar_result'
    freq_data_filename = os.path.join(
             phonts_sim_dir,
             'freq.dat')
    phon_lifetime_data_filename = os.path.join(
             phonts_sim_dir,
             'phon_lifetime.dat')

    bte_data = PhontsBteData(natoms=4,directory=phonts_sim_dir)
    bte_data.read()
    bte_data.build_data_at_temp(temp=400)
    ph_freq = bte_data.data[400][:,4]
    ph_lt = bte_data.data[400][:,5]
