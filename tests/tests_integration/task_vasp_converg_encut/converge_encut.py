import numpy as np
import pypospack.io.vasp as vasp
import pypospack.task.vasp as tsk_vasp
import os


if __name__ == '__main__':
    structure_filename = '../MgO_vasp_structuralminimization/MgO_calc/CONTCAR'
    xc = 'GGA'
    is_restart = True
    poscar = vasp.Poscar()
    poscar.read(structure_filename)

    # energy cutoff
    potcar = vasp.Potcar()
    potcar.symbols = poscar.symbols
    potcar.xc = xc
    potcar.write('POTCAR.tmp')
    potcar.read('POTCAR.tmp')
    os.remove('POTCAR.tmp')
    encut_min = max(potcar.encut_min)
    encut_max = 1.5 * max (potcar.encut_max)
    encut_step = 25 #eV

    print(encut_min,encut_max)

    incar_config = {\
        'ismear':0,
        'sigma':0.05,
        'ispin':1}

    slurm_dict = {\
            'email':'eragasa@ufl.edu',
            'qos':'phillpot',
            'ntasks':16,
            'time':"1:00:00"}
    tasks = {}
    task_results = {}

    for encut in np.arange(encut_min,encut_max+1,encut_step):
        encut = int(encut)
        str_encut = str(encut)
        tasks[str_encut] = tsk_vasp.VaspSimulation(
		    task_name = str_encut,
		    task_directory = str_encut,
		    restart=is_restart)

        incar_config['encut']=encut
        if tasks[str_encut].status == 'INIT':
            tasks[str_encut].config(poscar=structure_filename,incar=incar_config,xc=xc)

        if tasks[str_encut].status == 'CONFIG':
            tasks[str_encut].run(\
                job_type='slurm',
                exec_dict=slurm_dict)

        if tasks[str_encut].status == 'RUN':
            pass

        if tasks[str_encut].status == 'POST':
            tasks[str_encut].postprocess()
            total_energy = tasks[str_encut].outcar.total_energy
            n_atoms = len(tasks[str_encut].contcar.atomic_basis)
            total_energy_per_atom = total_energy / n_atoms

            task_results[str_encut] = {}
            task_results[str_encut]['encut'] = tasks[str_encut].outcar.encut
            task_results[str_encut]['total_energy'] = total_energy
            task_results[str_encut]['n_atoms'] = n_atoms
            task_results[str_encut]['total_energy_per_atom'] = total_energy / n_atoms

    for i,encut in enumerate(np.arange(encut_min,encut_max+1,encut_step)):
        encut = int(encut)
        str_encut = str(encut)
        str_old_encut = str(int(encut-25))
        if i == 0:
            print(str_encut,
                task_results[str_encut]['encut'],
                task_results[str_encut]['total_energy_per_atom']) 
        else:
            print(str_encut,
                task_results[str_encut]['encut'],
                task_results[str_encut]['total_energy_per_atom'],
                100*(task_results[str_encut]['total_energy_per_atom']-task_results[str_old_encut]['total_energy_per_atom']))
         

