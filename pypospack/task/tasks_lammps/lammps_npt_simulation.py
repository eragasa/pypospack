import os, random, sys
from collections import OrderedDict
from pypospack.task.lammps import LammpsSimulation
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal

class LammpsNptSimulation(LammpsSimulation):
    """ Class for LAMMPS structural minimization

    This data class defines additional attributes and methods necessary to 
    interact with the Workflow manager.

    Args:
        task_name(str): unique id for the task name being define
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS

    Attributes:
        config
        config_map
    """
    def __init__(self,
            task_name,
            task_directory,
            structure_filename,
            restart=None,
            fullauto=None,
            temperature=None,
            pressure=0,
            time_total=None,
            time_step=None,
            supercell=[10,10,10]):
       
        assert temperature is not None
        assert time_total is not None
        assert time_step is not None

        LammpsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                structure_filename=structure_filename)

        self.temperature = temperature
        self.pressure = pressure

        self.time_total = time_total
        self.time_step = time_step
        
        self.supercell = supercell

        self.lammps_out_fn = 'lammps.out'
        self.lattice_fn = 'lattice.out'
    
    def get_task_name(structure,temperature):
        task_name = '{s}.lmps_npt_{T}'.format(s=structure,T=str(int(T)))
        return task_name
    
    def postprocess(self):
        LammpsSimulation.postprocess(self)

    def on_post(self,configuration=None):
        self._get_lattice_parameter_from_lattice_out_file(
            filename=os.path.join(self.task_directory,self.lattice_fn),
            supercell=self.supercell)

        LammpsSimulation.on_post(self,configuration=configuration)

    def write_lammps_structure_filename(self):
        _supercell = self.supercell
        _structure_filename_src \
                = self.configuration['structure']['structure_filename']
        _structure = vasp.Poscar()
        _structure.read(_structure_filename_src)

        self.structure = crystal.make_super_cell(
                structure=_structure,
                sc=_supercell)

        self.lammps_structure = lammps.LammpsStructure(\
                obj=self.structure)

        _str_out = self.lammps_input_file_to_string()

        _lammps_filename = os.path.join(self.task_directory,filename)
        with open(_lammps_filename,'w') as f:
            f.write(_str_out)
    
    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_npt_thermostat(
                    time_total = self.time_total,
                    time_step = self.time_step,
                    temperature = self.temperature,
                    pressure=self.pressure),
                self._lammps_input_out_section(),
                ])
        return(str_out)

    def _lammps_input_create_atoms(self):
        str_out = LammpsSimulation._lammps_input_create_atoms(self)
        str_out += "replicate {} {} {}\n".format(
                self.supercell[0],
                self.supercell[1],
                self.supercell[2])
        str_out += "change_box all x scale 1 y scale 1 z scale 1 remap\n"

        return str_out

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'fix 1 all box/relax iso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe Lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            'unfix 1\n'
            )
        return str_out

    def _lammps_input_out_section(self):
        str_out = (
            '# ---- define output variables ----\n'
            'variable natoms equal "count(all)"\n'
            'variable tot_energy equal "c_eatoms"\n'
            'variable a11 equal "xhi-xlo"\n'
            'variable a22 equal "yhi-ylo"\n'
            'variable a33 equal "zhi-zlo"\n'
            'variable tilt_xy equal "xy"\n'
            'variable tilt_xz equal "xz"\n'
            'variable tilt_yz equal "yz"\n'
            'variable tot_press equal "press"\n'
            'variable press_xx equal "pxx"\n'
            'variable press_yy equal "pyy"\n'
            'variable press_zz equal "pzz"\n'
            'variable press_xy equal "pxy"\n'
            'variable press_xz equal "pxz"\n'
            'variable press_yz equal "pyz"\n'
            '\n'
            '# ---- output ----\n'
            'print \"pypospack:output_section:begin\"\n'
            'print \"num_atoms = ${natoms}"\n'
            'print \"a11 = ${a11}\"\n'
            'print \"a22 = ${a22}\"\n'
            'print \"a33 = ${a33}\"\n'
            'print \"a12 = ${tilt_xy}\"\n'
            'print \"a13 = ${tilt_xz}\"\n'
            'print \"a23 = ${tilt_yz}\"\n'
            'print \"tot_press = ${tot_press}\"\n'
            'print \"pxx = ${press_xx}\"\n'
            'print \"pyy = ${press_yy}\"\n'
            'print \"pzz = ${press_zz}\"\n'
            'print \"pxy = ${press_xy}\"\n'
            'print \"pxz = ${press_xz}\"\n'
            'print \"pyz = ${press_yz}\"\n'
            'print \"pypospack:output_section:done\"\n'
            'print \"pypospack:lammps_sim:done\"\n'
                  )
        return str_out
    def _determine_temperature_dampening(self,dt):
        # A Nose-Hoover thermostat will not work well for arbitrary values of 
        # Tdamp. If Tdamp is too small, the temperature can fluctuate wildly; 
        # if it is too large, the temperature will take a very long time to 
        # equilibrate. A good choice for many models is a Tdamp of around 100 
        # timesteps.
        # Ref: http://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command
        
        _tempdamp = dt*100
        return _tempdamp

    def _determine_pressure_dampening(self,dt):
        # A Nose-Hoover barostat will not work well for arbitrary values of 
        # Pdamp. If Pdamp is too small, the pressure and volume can fluctuate 
        # wildly; if it is too large, the pressure will take a very long time 
        # to equilibrate. A good choice for many models is a Pdamp of around 
        # 1000 timesteps. 
        # Ref: http://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command
        
        _pressdamp = dt*1000
        return _pressdamp
   
    def _determine_drag_coefficient(self):
        # In some cases (e.g. for solids) the pressure (volume) and/or 
        # temperature of the system can oscillate undesirably when a Nose/Hoover 
        # barostat and thermostat is applied. The optional drag keyword will 
        # damp these oscillations, although it alters the Nose/Hoover equations. 
        # A value of 0.0 (no drag) leaves the Nose/Hoover formalism unchanged. 
        # A non-zero value adds a drag term; the larger the value specified, 
        # the greater the damping effect. Performing a short run and monitoring 
        # the pressure and temperature is the best way to determine if the drag 
        # term is working. Typically a value between 0.2 to 2.0 is sufficient 
        # to damp oscillations after a few periods. Note that use of the drag 
        # keyword will interfere with energy conservation and will also change 
        # the distribution of positions and velocities so that they do not 
        # correspond to the nominal NVT, NPT, or NPH ensembles. 
        
        _drag = 1.0
        return _drag

    def _get_random_seed(self):

        _max_size = 10000000
        _seed = random.randrange(_max_size)
   
        return _seed

    def _lammps_input_npt_thermostat(self,
            time_step=0.001,
            time_total=100,
            temperature=100,
            pressure=0,
            seed1=None,
            seed2=None):
        """

        Args:
            time_step (float): time increment in picoseconds
            time_total (float): time simulation time in picoseconds
            temperature (int): simulation temperature in Kelvin
            pressure (int): simulation pressure in KBar
        """

        _dt = time_step
        _time1 =time_total
        _n_time_steps = int(_time1 / _dt)
        _temp0 = temperature
        _temp1 = temperature
        _press0 = pressure
        _press1 = pressure
        
        _tempdamp = self._determine_temperature_dampening(dt=_dt)
        _pressdamp = self._determine_pressure_dampening(dt=_dt)
        _drag = self._determine_drag_coefficient()
        
        if seed1 is None:
            _seed1 = self._get_random_seed()
        else:
            _seed1 = seed1

        if seed2 is None:
            _seed2 = self._get_random_seed()
        else:
            _seed2 = seed2

        str_out = "\n".join([
            "#------------------------------------------------------------------------------",
            "# RUN THERMOSTAT",
            "# running using an NPT Nose-Hoover style thermostat",
            "#------------------------------------------------------------------------------",
            "variable tempdamp equal {tempdamp}".format(tempdamp=_tempdamp),
            "variable pressdamp equal {pressdamp}".format(pressdamp=_pressdamp),
            "",
            "timestep {dt}".format(dt=_dt),
            "# set thermo -----------------------------------------------------------------",
            "thermo 100",
            "thermo_style custom step temp pe ke etotal press lx ly lz press pxx pyy pzz pxy pxz pyz vol",
            "thermo_modify flush yes",
            "# set averaging ---------------------------------------------------------------",
            "variable boxx equal lx",
            "variable boxy equal ly",
            "variable boxz equal lz",
            "variable boxp equal press",
            "variable boxt equal temp",
            "# calculate averages every 10000 steps",
            "# set initial velocities ------------------------------------------------------",
            "reset_timestep 0",
            "velocity all create {temp} {seed} mom yes dist gaussian loop all".format(
                temp=_temp0,
                seed=_seed1),
            "# fix for Nose-Hoover style thermostat ----------------------------------------",
            "fix 20 all npt temp {temp0} {temp1} {tempdamp} aniso 0.0 0.0 {pressdamp} drag {drag} couple xyz".format(
                temp0=_temp0,temp1=_temp1,tempdamp=_tempdamp,
                press0=_press0,press1=_press1,pressdamp=_pressdamp,
                drag=_drag),
            "fix 21 all ave/time 1 5000 5000 v_boxx v_boxy v_boxz v_boxp v_boxt file lattice.out",
            "run {n_time_steps}".format(n_time_steps=_n_time_steps),
            "unfix 20",
            "unfix 21",
        ])

        return str_out

    def _get_lattice_parameter_from_lattice_out_file(self,filename,supercell):
        with open(filename) as f:
            lines = f.readlines()

        names = [v.strip() for v in lines[1].strip().split(" ")[1:]]
        values = [float(v) for v in lines[len(lines)-1].strip().split(" ")]

        if self.results is None:
            self.results = OrderedDict()

        _task_name = self.task_name
        _a1 = values[names.index('v_boxx')]/supercell[0]
        _a2 = values[names.index('v_boxy')]/supercell[1]
        _a3 = values[names.index('v_boxz')]/supercell[2]

        self.results['{}.{}'.format(_task_name,'a1')] = _a1
        self.results['{}.{}'.format(_task_name,'a2')] = _a2
        self.results['{}.{}'.format(_task_name,'a3')] = _a3

        return self.results
