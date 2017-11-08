import os
from pypospack.task import Task
import pypospack.io.phonts as phonts
class PhontsSimulation(Task):
    """
    Args:
        task_name(str): This is the name of the simulation of the simulation.
            The primary purpose of this name is to be used as the jobname
            when submitting the job using a submission script, and also
            to store the unique key value of a simulation within a workflow.
        task_directory(str): This is the location of the directory for
            where the simulation files exist.
        structure_filename(str): This is the location of the POSCAR of the 
            structure which we want to study.
        restart(str): When set to True, this class will attempt to restart
            the PhonTS simulation from scratch
    Attributes:
        task_name(str): The name of the simulation job.
        task_directory(str): The path of the directory where this simulation
            is to be executed.
        structure_filename(str): Attribute initialized to None.
        structure(pypospack.crystal.SimulationCell): Attribute initialized to 
            None
        filename(str): name of the phonts input file  Initialized to be
            <task_dir>/phonons_input.dat

        potential(pypospack.potential.Potential)
        fp_interface(str): must be
            VASP, QE, or LAMMPS
        phonts_bin(str): the path of the phonts executable, by default 
            pypospack looks for PHONTS_BIN in the PATH os variable
        vasp_bin(str): the path of the vasp executable, by default
            pypoppack looks in VASP_BIN in the PATH of the os variable

    """
    def __init__(self,
            task_name='phonts',
            task_directory='phonts',
            structure_filename='POSCAR',
            phonts_filename='phonons_input.dat',
            restart=False):
        Task.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                restart=False)
        self.PHONT_INTERACTIONS = phonts.PHONTS_INTERACTIONS

        self.structure_filename = structure_filename
        self.structure = None
        self.phonts_filename = os.path.join(
                task_directory,phonts_filename)

        # external force evaluation
        self._interaction_type = None
        self._fp_interface = None
        self.is_prep_step = None
        self.numerical_2der = True
        self.numerical_3der = True
        self.d3_cutoff = None
        self.delta = 0.005
        self.phonts_bin = os.environ['PHONTS_BIN']

        # internal interatomic potential parameters
        # needed for passing empirical potentials --- to be implemented
        self.potential = None
        # job manager attributes
        self.job_scheduler = None
        self.slurm_phonts_dict = None
        self.slurm_vasp_dict = None

        # phonon calculation
        self.PHONTS_CALC_TYPES = [
                'phonons_only','iterations','qha']
        self._calculation_type = None
        self.iter_steps = 10
        self.mixing = 0.5
        self.qha = None 
        self.optimize = False
        
        self.phonon_kpoints = [6,6,6,10]
        self.pdos = [0.,30.,200,10.]

        # optimization keywords
        # phonts_calc_types
        self.quench_steps = 1000
        self.CELL_CONTROL_TYPES = [
                'fixed_cell', 'free_cell',
                'orthorhombic_cell', 'hydrostatic_pressure']
        self.cell_control = 'fixed_cell'
        self.atom_min_step = 0.000001
        self.cell_min_step = 0.1
        self.max_force = 1e-8
        self.quench_only = False
        self.dump = 1

        # set private attributes
        self._temperature = None
        # external conditions
        self.temperature = [300.,1200.,9]
        self.pressure = None
    @property
    def calculation_type(self):
        return self._calculation_type

    @calculation_type.setter
    def calculation_type(self,calculation_type):
         if calculation_type is 'phonons_only':
             self._calculation_type = calculation_type
         elif calculation_type is 'iterations':
             self._calculation_type = calculation_type
         elif calculation_type is 'qha':
             self._calcualtion_type = calculation_type
         else:
             msg_err = "calculation_type not supported"
             raise ValueError(msg_err)

    @property
    def interaction_type(self):
        return self._interaction_type

    @interaction_type.setter
    def interaction_type(self,interaction_type):
        if interaction_type == 'AbInitio':
            self._interaction_type = interaction_type
            self.fp_interface = 'VASP'
            self.numerical_2der = True
            self.numerical_3der = True
            self.delta = 0.005
        else:
            msg_err = "cannot set the PhonTS interface to {}"
            msg_err = msg_err.format(interaction_type)
            raise ValueError(msg_err)

    @property
    def fp_interface(self):
        return self._fp_interface

    @fp_interface.setter
    def fp_interface(self,fp_interface):
        if self.interaction_type != 'AbInitio':
            self.interaction_type = 'AbInitio'
        if fp_interface in phonts.PHONTS_FP_INTERFACES:
            self._fp_interface = fp_interface
        else:
            raise ValueError("invalid fp_interface parameter")

    @property 
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, temperature):
        """
        We can accept two forms of temperature:
        -- a numeric variable
        -- a list of numerics of length one --> [float(T)]
        -- a list of numerics of length three -->[float(t0),float(t1),int(N_t)]
        """
        if type(temperature) in [float,int]:
            self._temperature = [float(temperature)]
        elif type(temperature) is list:
            if len(temperature) == 1:
                self._temperature = list[temperature]
            elif len(temperature) == 3:
                self._temperature = [0,0,0]
                self._temperature[0] = float(temperature[0])
                self._temperature[1] = float(temperature[1])
                self._temperature[2] = int(temperature[2])
            else:
                msg_err = "temperature must be either a single temperature"
                msg_err += "or [float(t1),float(t2),int(N)]"
                raise ValueError(msg_err)
        else:
            msg_err = "temperature must be either a single temperature"
            msg_err += "or [float(t1),float(t2),int(N)]"
            raise ValueError(msg_err)

    def write_phonts_input_file(self,filename='phonons_input.dat'):
        raise NotImplementedError()
    
    def write_submission_script(self):

        if self.job_scheduler == 'slurm':
            self.slurm_phonts_dict['filename'] = os.path.join(
                    self.task_directory,
                    self.slurm_phonts_dict['filename'])
            slurm.write_phonts_batch_script(\
                    filename=self.slurm_phonts_dict['filename'],
                    job_name=self.slurm_phonts_dict['job_name'],
                    email=self.slurm_phonts_dict['email'],
                    qos=self.slurm_phonts_dict['qos'],
                    ntasks=self.slurm_phonts_dict['ntasks'],
                    time=self.slurm_phonts_dict['time'])
        elif self.job_scheduler == 'sge':
            raise NotImplementedError()
        elif self.job_scheduler == 'pbs':
            raise NotImplementedError()

    def run(self):
        self.orig_dir = os.getcwd()
        os.chdir(self.task_directory)
        if self.job_scheduler == 'slurm':
            cmd_s = 'sbatch {}'.format(self.slurm_phonts_dict['filename'])
            args = shlex.split(cmd_s)
            p = subprocess.Popen(args)

    def force_evaluation_to_string(self):
        str_out = ''
        if  self.fp_interface in phonts_fp_interfaces:
            str_out += self._get_external_force_evaluation_to_string() 
        elif self.potential is not None:
            str_out += self._get_internal_force_evaluation_to_string()
        else:
            err_msg = "either the potential must be specified or an "
            err_msg += "external method must be provided"
            raise ValueError(err_msg)
        return str_out
    
    def _get_external_force_evaluation_to_string(self):
        str_out = ''
        if self.fp_interface == 'VASP':
            pass
        elif self.fp_interface == 'QE':
            raise NotImplementedError('fp_interface == QE not supported')
        elif self.fp_interface == 'LAMMPS':
            raise NotImplementedError('fp_interface == LAMMPS not supported')
        return str_out

    def _get_internal_force_evaluation_to_string(self):
        str_out = ''
        str_out += ' '.join([str(v) for v in self.phonts_potential_type]) + '\n'
        str_out += ' '.join([str(v) for v in self.phonts_potential_params]) + '\n'
        return str_out 
    
    def simulation_cell_to_string(self):
        if isinstance(self.structure_filename,str):
            self.structure = vasp.Poscar()
            self.structure.read(self.structure_filename)

        str_out = ""
        if isinstance(self.structure,crystal.SimulationCell):
            # get relevant information
            n_species = len(self.structure.symbols)
            a0 = self.structure.a0
            a1 = self.structure.a1/a0
            a2 = self.structure.a2/a0
            a3 = self.structure.a3/a0
            n_atoms = self.structure.n_atoms

            # build string
            str_out += "species {}\n".format(n_species)
            for s in self.structure.symbols:
                s_charge = self._get_charge(s)
                s_amu = crystal.get_amu(s)
                str_out += "{} {} {}\n".format(s,s_amu,s_charge)
            str_out += "Lattice {:.10f}\n".format(a0)
            str_out += "cell {:.10f} {:.10f} {:.10f}\n".format(a1,a2,a3)
            str_out += "natoms {}\n".format(n_atoms)
            str_out += "fractional\n"
            for a in self.structure.atomic_basis:
                s = a.symbol
                x = a.position[0]
                y = a.position[1]
                z = a.position[2]
                str_out += "{} {:10f} {:10f} {:10f}\n".format(s,x,y,z)
        else:
            msg_out = "structure must be an instance of "
            msg_out += "pypospack.crystal.SimulationCell"
            raise ValueError(msg_out)

        return str_out

class PhontsBte(PhontsSimulation):
    def __init__(self,
            task_name,
            task_directory,
            phonts_filename,
            structure_filename,
            restart=True):
        PhontsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                phonts_filename=phonts_filename,
                structure_filename=structure_filename,
                restart=restart)
        
        self._interaction_type='AbInitio'
        self.fp_interface='VASP'
        self.numerical_2der=True
        self.numerical_3der=True
        self.d3_cutoff=5.0
        self.delta=0.005

        self.calculation_type = 'iterations'
        self.iter_steps = 10
        self.phonon_kpoints = [6,6,6,10]
        self.temperature = [300, 1200,9]
    
    def write_phonts_input_file(self):
        phonts_config_dict = OrderedDict()
        phonts_config_dict['interaction_type'] = self.interaction_type
        phonts_config_dict['is_prep_step'] = self.is_prep_step
        phonts_config_dict['fp_interface'] = self.fp_interface
        phonts_config_dict['numerical_2der'] = self.numerical_2der
        phonts_config_dict['numerical_3der'] = self.numerical_3der
        phonts_config_dict['d3_cutoff'] = self.d3_cutoff
        phonts_config_dict['delta'] = self.delta
     
        phonts_config_dict['calculation_type'] = self.calculation_type
        phonts_config_dict['iter_steps'] = self.iter_steps
        phonts_config_dict['phonon_kpoints'] = self.phonon_kpoints
        phonts_config_dict['temperature']= self.temperature
        phonts_input_file = PhontsInputfile(
                filename=self.phonts_filename,
                structure_filename=self.structure_filename,
                phonts_config_dict=phonts_config_dict)

class PhontsBte_PrepStep(PhontsSimulation):
    def __init__(self,
            task_name,
            task_directory,
            simulation_filename,
            phonts_filename,
            restart=True):
        PhontsBte.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                phonts_filename=phonts_filename,
                structure_filename=structure_filename,
                restart=restart)
        
        self.is_prep_step = True

class PhontsBte_IntermediateStep(PhontsSimulation):
    def __init__(self,
            task_name,
            task_directory,
            simulation_filename,
            phonts_filename,
            restart=True):
        PhontsBte.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                phonts_filename=phonts_filename,
                structure_filename=structure_filename,
                restart=restart)

class PhontsBte_FinalStep(PhontsSimulation):
    def __init__(self,
            task_name,
            task_directory,
            simulation_filename,
            phonts_filename,
            restart=True):
        PhontsBte.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                phonts_filename=phonts_filename,
                structure_filename=structure_filename,
                restart=restart)
      
        self.is_prep_step = False
