import random, sys
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
            structure_filename):
        LammpsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                structure_filename=structure_filename)

    def postprocess(self):
        LammpsSimulation.postprocess(self)

    def write_lammps_structure_filename(self):
        _supercell = self.configuration['structure']['supercell']
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
        _dt = self.configuration['task']['dt']
        _T = self.configuration['task']['temperature']
        
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section(),
                self._lammps_input_nvt_thermostat(dt=_dt,temp=_T)
                ])
        return(str_out)

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
            # 'thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            'unfix 1\n'
            )
        return str_out

    def _lammps_input_nvt_thermostat(self,dt=0.001,temp=100,press=0,seed=None):
        _dt = dt
        _time1 = 100
        _n_time_steps = int(_time1 / _dt)
        _temp0 = temp
        _temp1 = temp
        _press0 = press
        _press1 = press
        # A Nose-Hoover thermostat will not work well for arbitrary values of 
        # Tdamp. If Tdamp is too small, the temperature can fluctuate wildly; 
        # if it is too large, the temperature will take a very long time to 
        # equilibrate. A good choice for many models is a Tdamp of around 100 
        # timesteps.
        # Ref: http://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command
        _tempdamp = dt*100
        # A Nose-Hoover barostat will not work well for arbitrary values of 
        # Pdamp. If Pdamp is too small, the pressure and volume can fluctuate 
        # wildly; if it is too large, the pressure will take a very long time 
        # to equilibrate. A good choice for many models is a Pdamp of around 
        # 1000 timesteps. 
        _pressdamp = dt*1000
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
        if seed is None:
            # this doesn't work because python maxsize is larger than C int
            #_max_size = sys.maxsize
            _max_size = 10000000
            _seed = random.randrange(_max_size)
          
        else:
            _seed = seed
        str_out = "\n".join([
"#------------------------------------------------------------------------------",
"# RUN THERMOSTAT",
"# running using an NPT Nose-Hoover style thermostat",
"#------------------------------------------------------------------------------",
"variable tempdamp equal {tempdamp}".format(tempdamp=_tempdamp),
"variable pressdamp equal {pressdamp}".format(pressdamp=_pressdamp),
"# set initial velocities",
"velocity all create {temp} {seed} mom yes dist gaussian loop all".format(
    temp=_temp0,seed=_seed),
"# fix for Nose-Hoover style thermostat",
"fix 20 all npt temp {temp0} {temp1} {tempdamp} aniso 0.0 0.0 {pressdamp} drag {drag} couple xyz".format(
    temp0=_temp0,temp1=_temp1,tempdamp=_tempdamp,
    press0=_press0,press1=_press1,pressdamp=_pressdamp,
    drag=_drag),
"variable boxx equal lx",
"variable boxy equal ly",
"variable boxz equal lz",
"variable boxp equal press",
"variable boxt equal temp",
"# calculate averages every 10000 steps",
"reset_timestep 0",
"fix 21 all ave/time 10 10 100 v_boxx v_boxy v_boxz v_boxp v_boxt file lattice.out",
"thermo_style custom step temp pe ke etotal press lx ly lz press pxx pyy pzz pxy pxz pyz vol",
"thermo_modify flush yes",
"thermo 100",
"run {n_time_steps}".format(n_time_steps=_n_time_steps),
"unfix 20",
"unfix 21",
"print \"seeds:{seed}\"".format(seed=_seed)])

        return str_out

