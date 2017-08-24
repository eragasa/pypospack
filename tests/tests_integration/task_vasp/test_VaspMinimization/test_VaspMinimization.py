import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

def vasp_structural_minimization(structure,task_name):
    assert isinstance(structure,pypospack.crystal.SimulationCell)
    assert isinstance(task_name,str)

    vasp_simulation = vasp.VaspSimulation()
    vasp_simulation.simulation_directory = task_dir
    vasp_simulation.poscar = vasp.Poscar(structure)
    vasp_simulation.incar = vasp.Incar()
    vasp_simulation.potcar = vasp.Potcar()
    vasp_simulation.kpoints = vasp.Kpoints()
def stuff():
    print('create_structure')


if __name__ == "__main__":
    symbol = 'Ni'
    sg = 'fcc'
    a = 3.508
    task_name = 'minimize_init'
    job_name = "{}_{}_min1".format(symbol,sg)
    email= "eragasa@ufl.edu"
    qos="phillpot"
    ntasks=16
    time="1:00:00"

    init_directory = os.getcwd()
    task_dir = os.path.join(os.getcwd(),task_name)

    vasp_simulation = vasp.VaspSimulation()
    vasp_simulation.simulation_directory = task_dir
    vasp_simulation.poscar = vasp.Poscar(\
            ase.build.bulk(symbol,sg,a,cubic=True))
    vasp_simulation.incar = vasp.Incar()
    vasp_simulation.potcar = vasp.Potcar()
    vasp_simulation.kpoints = vasp.Kpoints()

    # create directory
    # create vasp simulation directory and files
    # delete the directory and the contents if it exists
    if os.path.exists(vasp_simulation.simulation_directory):
        shutil.rmtree(vasp_simulation.simulation_directory)
    os.mkdir(vasp_simulation.simulation_directory)

    # configure potcar file
    print('confguring potcar file...')
    vasp_simulation.xc = 'GGA'
    vasp_simulation.symbols = vasp_simulation.poscar.symbols
    vasp_simulation.potcar.symbols = vasp_simulation.symbols
    vasp_simulation.potcar.xc = vasp_simulation.xc

    # get information from potcar
    vasp_simulation.potcar.write(
            os.path.join(vasp_simulation.simulation_directory,"POTCAR"))
    vasp_simulation.potcar.read(
            os.path.join(vasp_simulation.simulation_directory,"POTCAR"))

    # stuff to put into tests
    # use the largest enmax of all a
    # print('symbols:{}'.format(','.join(vasp_simulation.potcar.symbols)))
    # print('recommended encut range: {}-{}'.format(
    #     min(vasp_simulation.potcar.encut_min)
    #     max(vasp_simulation.potcar.encut_max)
    #    vasp_simulation.encut_min,
    #    vasp_simulation.encut_max))
    # print('n_atoms:{}'.format(vasp_simulation.poscar.n_atoms))
    # print('n_atomic_basis:{}'.format(len(vasp_simulation.poscar.atomic_basis)))
    # print('n_interstitials:{}'.format(len(vasp_simulation.poscar.interstitials)))
    # configure incar file
    
    vasp_simulation.encut = max(vasp_simulation.potcar.encut_max)
    vasp_simulation.natoms = vasp_simulation.poscar.n_atoms
    vasp_simulation.incar.encut = vasp_simulation.encut

    magmom_init = 1.0
    vasp_simulation.incar.ismear = 1      # methfessel paxton order 1
    vasp_simulation.incar.sigma = 0.20    # smearing parameter

    vasp_simulation.incar.ispin = 2       # spin-polarized calculations
    vasp_simulation.incar.magmom = "{}*{}".format(
            vasp_simulation.poscar.n_atoms,magmom_init)
    vasp_simulation.incar.ibrion = 2      # conjugate gradient method
    vasp_simulation.incar.isif = 3        # relax everything
    vasp_simulation.incar.potim = 0.5     # dampening parameter
    vasp_simulation.incar.ediffg = -0.001 # ev/Ang

    vasp_simulation.poscar.write(
            os.path.join(vasp_simulation.simulation_directory,"POSCAR"))
    vasp_simulation.incar.write(
            os.path.join(vasp_simulation.simulation_directory,"INCAR"))
    vasp_simulation.kpoints.write(
            os.path.join(vasp_simulation.simulation_directory,"KPOINTS"))

    slurm.write_vasp_batch_script(
            filename=os.path.join(
                vasp_simulation.simulation_directory,"runjob_hpg.slurm"),
            job_name=job_name,
            email=email,
            qos=qos,
            ntasks=ntasks,
            time=time)
    
    init_dir = os.getcwd()
    os.chdir(vasp_simulation.simulation_directory)
    result = subprocess(sbatch runjob_hpg.slurm)
    with open(job.info,'w') as fout:
        fout.write(result.stdout)
    fout.close()
