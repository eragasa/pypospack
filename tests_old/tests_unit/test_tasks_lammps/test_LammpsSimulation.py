import pypospack.task.lammps as tsk_lammps
import pypospack.potential as potential

if __name__ == '__main__':
    task_name = 'task_name'
    task_directory = 'task_name'
    task = tsk_lammps.LammpsSimulation(\
            task_name = task_name,
            task_directory = task_directory)

    print('task.status:{}'.format(task.status))

    pot = potential.Buckingham(['Mg','O'])
    task.potential = pot
    print('write_lammps_input_file():{}'.format(task.lammps_input_file_to_string()))
