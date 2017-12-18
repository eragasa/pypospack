import os
import numpy as np

class LammpsNptSimulationData(object):

    def __init__(self,
            task_directory):
        self.task_directory = task_directory
        self.latticeout_filename = os.path.join(
                self.task_directory,
                'lattice.out')

        self.lattice_labels = None
        self.lattice_data = None

    def read_lattice_out(self):
        with open(self.latticeout_filename,'r') as f:
            lines = f.readlines()

        _lattice_labels = None
        _lattice_data = []
        for i,line in enumerate(lines):
            if i in [0]:
                pass
            elif i in [1]:
                _args = line.split('#')[1].split(' ')
                _lattice_labels = [a.strip() for a in _args]
                _lattice_labels = [l for l in _lattice_labels if l != '']
            else:
                _args = line.split(' ')
                _args = [float(a) for a in _args]
                _lattice_data.append(list(_args))

        self.lattice_labels = list(_lattice_labels)
        self.lattice_data = np.array(_lattice_data)
        self.lattice_data[4] = self.lattice_data[4]/10000 # bar to GPa
    def plot_lattice_data(self):
        if self.lattice_labels is None or self.lattice_data is None:
            self.read_lattice_out()

        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3,1)
        axes[0].plot(self.lattice_data[:,0],self.lattice_data[:,1])
        axes[1].plot(self.lattice_data[:,0],self.lattice_data[:,4])
        axes[2].plot(self.lattice_data[:,0],self.lattice_data[:,5])
        #axes.plot(self.lattice_data[0],self.lattice_data[2])
        #axes.plot(self.lattice_data[0],self.lattice_data[3])
        fig.savefig("{}.png".format(self.task_directory))


task_directory = 'MgO_NaCl_unit.npt.T1000'
lammps_results = LammpsNptSimulationData(
        task_directory=task_directory)
lammps_results.read_lattice_out()
lammps_results.plot_lattice_data()
print(lammps_results.lattice_data[:,0])
