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

        for i,T in enumerate(self.temperatures):

            _task_name = "{}.npt.T{}".format(structure_prefix,T)
            _task_directory = os.path.join(
                    self.workflow_directory,
                    _task_name)
            self.lammps_npt_data[T] = LammpsNptSimulationData(
                    task_directory=_task_directory)
            self.lammps_npt_data[T].read_lattice_out()

    def calculate_thermal_expansion_curve(self,N=100):
        _labels = ['T','lx','ly','lz','P','T']
        _lammps_out = "lammps.out" 
        _th_exp_curve = []
        
        # get 0K
        _task_directory = "{}.npt.T{}".format(
                self.structure_prefix,
                self.temperatures[0])
        _lammps_out_filename = os.path.join(
                self.workflow_directory,
                _task_directory,
                _lammps_out)
        _lines = None
        with open(_lammps_out_filename,'r')as f:
            _lines = f.readlines()
        _row = [0,None,None,None,None,0]
        for _line in _lines:
            if _line.startswith('xx = '):
                _row[1] = float(_line.split('=')[1].strip())
            elif _line.startswith('yy = '):
                _row[2] = float(_line.split('=')[1].strip())
            elif _line.startswith('zz = '):
                _row[3] = float(_line.split('=')[1].strip())
            elif _line.startswith('tot_press ='):
                _row[4] = float(_line.split('=')[1].strip())
            elif all([v is not None for v in _row]):
                _th_exp_curve.append(_row)
                break
            else:
                pass

        for T in self.temperatures:
            _task_directory = "{}.npt.T{}".format(
                    self.structure_prefix,T)
            self.lammps_npt_data[T].get_simple_moving_averages(N=N)
            _last_sma = [s[-1] for s in self.lammps_npt_data[T].sma[N]]
            _th_exp_curve.append([T] + _last_sma)

        _fname_out = os.path.join(
                self.workflow_directory,
                'th_exp_curve.out')
        
        _temperatures = [0] + self.temperatures
        with open(_fname_out,'w') as f:
            f.write(" ".join(_labels) + "\n")
            for T in _temperatures:
                _line = ' '.join([str(row) for row in _th_exp_curve])
                f.write(_line + "\n")

        self.thermal_expansion_curve = np.array(_th_exp_curve)
        return _th_exp_curve

from pypospack.pyposmat import PyposmatDataFile

pyposmat_data_filename = os.path.join(
        'test_WorkflowLammpsThermalExpansion',
        'subselect.d_metric.sum_b_lt_median.out')
datafile = PyposmatDataFile(filename=pyposmat_data_filename)
datafile.read()
best_median_id = int(datafile.df.loc[datafile.df['sum_b_lt_median'].idxmax(axis=1),'sim_id'])
best_median_name = 'param_{}'.format(best_median_id)
best_dmetric_id = int(datafile.df.loc[datafile.df['d_metric'].idxmin(axis=1),'sim_id'])
best_dmetric_name = 'param_{}'.format(best_dmetric_id)
workflow_directories = ['MgO_LC','MgO_BG1','MgO_BG2']
workflow_directories += [v for v in os.listdir() if v.startswith('param_')]
workflow_data = OrderedDict()
for wfd in workflow_directories:
    workflow_data[wfd] = WorkflowLammpsThermalExpansionData(workflow_directory=wfd)
    workflow_data[wfd].read_data()
    workflow_data[wfd].calculate_thermal_expansion_curve(N=200)
    print(wfd,type(workflow_data[wfd]))


import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
workflow_color = OrderedDict()
for wfd in workflow_directories:
    workflow_color[wfd] = 'lightgrey'
workflow_color['MgO_LC'] = 'red'
workflow_color['MgO_BG1'] = 'purple'
workflow_color['MgO_BG2'] = 'blue'
workflow_color[best_median_name] = 'darkgrey'
workflow_color[best_dmetric_name] = 'black'

workflow_linewidth = OrderedDict()
for wfd in workflow_directories:
    workflow_linewidth[wfd] = 2.0
workflow_linewidth['MgO_LC'] = 4.0
workflow_linewidth['MgO_BG1'] = 4.0
workflow_linewidth['MgO_BG2'] = 4.0
workflow_linewidth[best_median_name] = 4.0
workflow_linewidth[best_dmetric_name] = 4.0

workflow_label = OrderedDict()
workflow_label['MgO_LC'] = 'LC+2.0'
workflow_label['MgO_BG1'] = 'BG+2.0'
workflow_label['MgO_BG2'] = 'BG+1.7'
workflow_label[best_median_name] = 'median'
workflow_label[best_dmetric_name] = 'dmetric'

_xlabel_str = r'Temperature [K]'
_ylabel_str = r'Lattice Parameter [\AA]'
fig,ax = plt.subplots(1,1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for i,wfd in enumerate(workflow_directories):
    print(wfd)
    if wfd in workflow_label:
        ax.plot(
                workflow_data[wfd].thermal_expansion_curve[:,0],
                workflow_data[wfd].thermal_expansion_curve[:,1],
                label=workflow_label[wfd],
                linewidth=workflow_linewidth[wfd],
                color=mcolors.cnames[workflow_color[wfd]])
    else:
        ax.plot(
            workflow_data[wfd].thermal_expansion_curve[:,0],
            workflow_data[wfd].thermal_expansion_curve[:,1],
            linewidth=workflow_linewidth[wfd],
            color=mcolors.cnames[workflow_color[wfd]])
ax.set_xlabel(_xlabel_str)
ax.set_ylabel(_ylabel_str)
legend = plt.legend()
fig.savefig('th_exp_curve_compare.png')
