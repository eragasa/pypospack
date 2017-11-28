import os,copy
from collections import OrderedDict
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
        self.sma = None

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
        #axes.plot(self.lattice_data[0],self.lattice_data[2])
        #axes.plot(self.lattice_data[0],self.lattice_data[3])
        axes[1].plot(self.lattice_data[:,0],self.lattice_data[:,4])
        axes[2].plot(self.lattice_data[:,0],self.lattice_data[:,5])
        fig.savefig("{}.png".format(self.task_directory))

    def get_simple_moving_averages(self,N=30):
        _kernel = np.ones((N,))/N
        _sma = [
                np.convolve(self.lattice_data[:,1],_kernel, mode='valid'),
                np.convolve(self.lattice_data[:,2],_kernel, mode='valid'),
                np.convolve(self.lattice_data[:,3],_kernel, mode='valid'),
                np.convolve(self.lattice_data[:,4],_kernel, mode='valid'),
                np.convolve(self.lattice_data[:,5],_kernel, mode='valid')]
        if self.sma is None:
            self.sma = OrderedDict()
            self.sma[N] = copy.deepcopy(_sma)
if False:
    T = 1000
    task_directory = 'MgO_LC/MgO_NaCl_unit.npt.T1000'
    lammps_results = LammpsNptSimulationData(
            task_directory=task_directory)
    lammps_results.read_lattice_out()
    lammps_results.plot_lattice_data()
    sma = lammps_results.get_simple_moving_averages(N=100)
    last_sma = [s[-1] for s in sma]
    print(' '.join([str(v) for v in ([T] + last_sma)]))

class WorkflowLammpsThermalExpansionData(object):

    def __init__(self,
            workflow_directory):

        self.workflow_directory = workflow_directory
        self.task_directories = None
        self.lammps_npt_data = None
    def read_data(self,structure_prefix='MgO_NaCl_unit',
            temperature_low=100,
            temperature_high=2000,
            temperature_step=100):

        self.structure_prefix = 'MgO_NaCl_unit'
        _tlow = int(temperature_low)
        _thigh = int(temperature_high)
        _tstep = int(temperature_step)

        self.temperatures = [t for t in range(_tlow,_thigh+1,_tstep)]
        
        self.task_directories = []
        
        self.lammps_npt_data = OrderedDict()
        for T in self.temperatures:
            _task_directory = "{}.npt.T{}".format(
                    structure_prefix,T)
            self.lammps_npt_data[T] = LammpsNptSimulationData(
                    task_directory=_task_directory)
            self.lammps_npt_data[T].read_lattice_out()

    def calculate_thermal_expansion_curve(self,N=100):
        _labels = ['T','lx','ly','lz','P','T']
        for T in self.temperatures:
            _task_directory = "{}.npt.T{}".format(
                    self.structure_prefix,T)
            self.lammps_npt_data[T].get_simple_moving_averages(N=N)

        _fname_out = os.path.join(
                self.workflow_directory,
                'th_exp_curve.out')
        
        with open(_fname_out,'w') as f:
            f.write(" ".join(_labels) + "\n")
            for T in self.temperatures:
                _last_sma = [s[-1] for s in self.lammps_npt_data[T].sma[N]]
                _line = ' '.join([str(v) for v in ([T] + _last_sma)])
                print(_line)
                f.write(_line + "\n")

_workflow_directory = 'MgO_LC'
workflow_data = WorkflowLammpsThermalExpansionData(
        workflow_directory='MgO_LC')
workflow_data.read_data(
        )
workflow_data.calculate_thermal_expansion_curve(N=200)

