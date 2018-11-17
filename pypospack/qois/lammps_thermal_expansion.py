from collections import OrderedDict
import numpy as np
from pypospack.qoi import Qoi

class ThermalExpansion(Qoi):

    """

    Args:
        temperature_min (float,int): beginning of the temperature range in Kelvin
        temperature_max (float,int): end of the temperature range in Kelvin
        temperature_step (float,int): increments of the temperature range in Kelvin
        time_total (int): total simulation time in fs
        time_step (int): simulation time step in fs
    """
    def __init__(self,qoi_name,structures,
            temperature_min=0,
            temperature_max=2700,
            temperature_step=100,
            time_total=10,
            time_step=0.001,
            supercell=[5,5,5]):

        _qoi_name = qoi_name
        _qoi_type = 'lmps_thermal_expansion'

        _structures = OrderedDict()
        _structures['ideal'] = structures['ideal']

        Qoi.__init__(self,
                qoi_name=_qoi_name,
                qoi_type=_qoi_type,
                structures=_structures)
    
        self.temperature_min = temperature_min
        self.temperature_max = temperature_max
        self.temperature_step = temperature_step
   
        self.time_total=time_total
        self.time_step=time_step

        self.supercell = supercell
    
    def determine_tasks(self):

        T = self.temperature_min
        
        while T <= self.temperature_max:
            if T == 0:
                _ideal_structure_name = self.structures['ideal']
                _ideal_task_type = 'lmps_min_all'
                _ideal_task_name = '{}.{}'.format(
                        _ideal_structure_name,
                        _ideal_task_type
                        )
                _bulk_structure_name = None

                self.add_task(
                        task_type=_ideal_task_type,
                        task_name=_ideal_task_name,
                        task_structure=_ideal_structure_name,
                        bulk_structure_name=_bulk_structure_name,
                        )
            else:
                _ideal_structure_name = self.structures['ideal']
                _ideal_task_type = 'lmps_npt'.format(T)
                _ideal_task_name = '{}.{}_{}'.format(
                        _ideal_structure_name,
                        _ideal_task_type,
                        T
                        )
                _bulk_structure_name = None

                self.add_task(
                        task_type=_ideal_task_type,
                        task_name=_ideal_task_name,
                        task_structure=_ideal_structure_name,
                        bulk_structure_name=_bulk_structure_name,
                        task_options={
                            'temperature':T,
                            'time_total':self.time_total,
                            'time_step':self.time_step,
                            'supercell':self.supercell}
                        )

            T = T + self.temperature_step

    def calculate_thermal_expansion_coefficient(self,temperatures,lattice_constants):

        assert isinstance(temperatures,list)
        assert isinstance(lattice_constants,list)
        T = list(temperatures)
        a0 = list(lattice_constants)

        a0_at_0K = a0[0]
        for i,v in enumerate(a0):
            a0[i] = v/a0_at_0K-1

        
        T = np.array(temperatures)
        a0 = np.array(a0)

        print(T)
        print(a0)

        T = T[:,np.newaxis]            # T needs to be a column vector
        # model is y = a*x
        alpha_L,_,_,_ = np.linalg.lstsq(T,a0)
        print('alpha_L:{}'.format(alpha_L[0]))
        return alpha_L[0]
    
    def calculate_qois(self,task_results):
        
        _prefix = '{}.{}'.format(self.structures['ideal'],self.qoi_type)
        s = self.structures['ideal']
        T = self.temperature_min

        lattice_constants = OrderedDict()
        while T <= self.temperature_max:
            if T == 0:
                lattice_constants[T] = np.sqrt(
                    task_results["{}.lmps_min_all.a11".format(s)]**2 \
                    + task_results['{}.lmps_min_all.a12'.format(s)]**2 \
                    + task_results["{}.lmps_min_all.a13".format(s)]**2
                )
            else:

                try:
                    lattice_constants[T] = task_results["{}.lmps_npt_{}.a1".format(s,T)]
                except KeyError as e:
                    for k,v in task_results.items():
                        print(k,v)
                    raise

            T = T + self.temperature_step

        self.qois = OrderedDict()
        
        # add lattice constants at different temperatures
        for k,v in lattice_constants.items():
            self.qois['{}.a0_{}'.format(_prefix,T)] = v
       
        _temperatures = [k for k,v in lattice_constants.items()]
        _lattice_constants = [v for k,v in lattice_constants.items()]
        
        # add thermal expansion coefficient
        self.qois['{}.thermal_expansion_coefficient'.format(_prefix)] = \
            self.calculate_thermal_expansion_coefficient(
                    temperatures=_temperatures,
                    lattice_constants=_lattice_constants
            )
