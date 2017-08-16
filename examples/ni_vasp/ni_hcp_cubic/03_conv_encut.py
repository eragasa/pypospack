import os, shutil
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm


if __name__ == "__main__":
    # DiRECTORY
    dir_min_0 = 'minimize_init'
    dir_conv_kpoints = 'conv_kpoints'
    dir_conv_encut = 'conv_encut'
    dir_min_f = 'minimize_f'
    
    # SLURM INFO
    slurm_email = "eragasa@ufl.edu"
    slurm_qos = "phillpot"
    slurm_ntasks = 16
    slurm_time = "1:00:00"
    
    # initial structure info
    symbol = 'Ni'
    sg = 'fcc'
    a = 3.508
    
    # potcar config
    xc = 'GGA'
    
    # INITIAL MINIMIZATION
    
    # get information from previous simulation

    # create directory to hold kpoint convergence simulations
    if os.path.exists(dir_conv_encut):
        shutil.rmtree(dir_conv_encut)
    os.mkdir(dir_conv_encut)

    kpcv_kp_mesh_size = [8,8,8]
    kpcv_encut_low = 250
    kpcv_encut_high = 600
    kpcv_encut_step = 25
    # do kpoint mesh
    for encut in range(kpcv_encut_low,kpcv_encut_high+1,kpcv_encut_step):
        vasp_simulation = vasp.VaspSimulation()
        vasp_simulation.simulation_directory = os.path.join(
                dir_conv_encut,
                '{}eV'.format(encut))
        slurm_job_name = "{}_{}_{}eV".format(
            symbol,sg,encut)

        if os.path.exists(vasp_simulation.simulation_directory):
            shutil.rmtree(vasp_simulation.simulation_directory)
        os.mkdir(vasp_simulation.simulation_directory)

        # write poscar
        vasp_simulation.poscar = vasp.Poscar(vasp.read_poscar_file(\
                filename = os.path.join(dir_min_0,'CONTCAR')))
        vasp_simulation.poscar.write(
                os.path.join(vasp_simulation.simulation_directory,"POSCAR"))

        # write potcar
        vasp_simulation.xc = xc
        vasp_simulation.symbols = vasp_simulation.poscar.symbols
        vasp_simulation.potcar.symbols = vasp_simulation.symbols
        vasp_simulation.potcar.xc = vasp_simulation.xc
        vasp_simulation.potcar.write(
                os.path.join(vasp_simulation.simulation_directory,"POTCAR"))

        # write incar
        vasp_simulation.incar = vasp.Incar()
        vasp_simulation.incar.read(\
                os.path.join(dir_min_0,'INCAR'))
        # change encut
        vasp_simulation.incar.encut = encut
        # turn off ionic relaxation
        vasp_simulation.incar.ibrion = None
        vasp_simulation.incar.isif = None
        vasp_simulation.incar.write(\
                os.path.join(
                    vasp_simulation.simulation_directory,
                    'INCAR'))
        # write kpoints
        vasp_simulation.kpoints = vasp.Kpoints()
        vasp_simulation.kpoints.meshsize= kpcv_kp_mesh_size
        vasp_simulation.kpoints.write(\
                os.path.join(
                    vasp_simulation.simulation_directory,
                    'KPOINTS'))

        # write slurm script
        slurm.write_vasp_batch_script(
            filename=os.path.join(\
                    vasp_simulation.simulation_directory,
                    "runjob_hpg.slurm"),
            job_name=slurm_job_name,
            email=slurm_email,
            qos=slurm_qos,
            ntasks=slurm_ntasks,
            time=slurm_time)

        # submit job 	
        init_dir = os.getcwd()
        os.chdir(vasp_simulation.simulation_directory)
        print('running in {}'.format(os.getcwd()))
        os.system('sbatch runjob_hpg.slurm')
        os.chdir(init_dir)
