import numpy as np
import pypospack.io.vasp as vasp
import pypospack.task.vasp as tsk_vasp
import os
if __name__ == '__main__':
    structure_filename = 'rsrc/MgO_NaCl_prim.vasp'
    xc = 'GGA'
    encut =  800
    poscar = vasp.Poscar()
    poscar.read(structure_filename)

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
    kpoint_meshes = [ [4,4,4],[6,6,6],[8,8,8],[10,10,10],[11,11,11],[12,12,12],[13,13,13]]
    for mesh_size in kpoint_meshes:
        task_name = '{}x{}x{}'.format(mesh_size[0],mesh_size[1],mesh_size[2])
        tasks[task_name]= tsk_vasp.VaspSimulation(
		    task_name = task_name,
		    task_directory = task_name,
		    restart=True)

        incar_config['encut']=encut
        kpoints_config = {'mesh_size':mesh_size}
        if tasks[task_name].status == 'INIT':
            tasks[task_name].config(
                poscar=structure_filename,
                incar=incar_config,
                kpoints=kpoints_config,xc=xc)

        if tasks[task_name].status == 'CONFIG':
            tasks[task_name].run(\
                job_type='slurm',
                exec_dict=slurm_dict)

        if tasks[task_name].status == 'RUN':
            pass

        if tasks[task_name].status == 'POST':
            tasks[task_name].postprocess()
            total_energy = tasks[task_name].outcar.total_energy
            n_atoms = len(tasks[task_name].contcar.atomic_basis)
            total_energy_per_atom = total_energy / n_atoms
            task_results[task_name] = {}
            task_results[task_name]['encut'] = tasks[task_name].outcar.encut
            task_results[task_name]['total_energy'] = total_energy
            task_results[task_name]['n_atoms'] = n_atoms
            task_results[task_name]['total_energy_per_atom'] = total_energy / n_atoms
    for i,mesh_size in enumerate(kpoint_meshes):
        if i == 0:
            t1 = '{}x{}x{}'.format(mesh_size[0],mesh_size[1],mesh_size[2])
            energy_per_atom = task_results[t1]['total_energy_per_atom']
            print(t1,energy_per_atom)
        else:
            t1 = '{}x{}x{}'.format(mesh_size[0],mesh_size[1],mesh_size[2])
            t0 = '{}x{}x{}'.format(kpoint_meshes[i-1][0],kpoint_meshes[i-1][1],kpoint_meshes[i-1][2])
            energy_per_atom = task_results[t1]['total_energy_per_atom']
            d_energy_per_atom = 1000*(task_results[t1]['total_energy_per_atom'] - task_results[t0]['total_energy_per_atom'])
            print(t1,energy_per_atom,d_energy_per_atom) 

