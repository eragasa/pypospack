import os, copy, shutil, subprocess
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

class GulpPhononCalculation(object):
    def __init__(self,taskname,task_dir):
        self.taskname = taskname
        self.task_directory = task_dir
        self.structure_filename = None
        self.simulation_cell = vasp.Poscar()
        self.filename = os.path.join(\
                self.task_directory,
                'gulp.in')
        self.potential = None

        self._create_task_directory()

    def _create_task_directory(self):
        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)
        os.mkdir(self.task_directory)

    def run(self):
        # write gulp input file
        try:
            self.write_gulp_input_file(\
                    filename=os.path.join(self.task_directory,'gulp.in'),
                    structure_filename=os.path.join(\
                            self.task_directory,
                            self.structure_filename))
        except:
            print("filename:{}".format(
                os.path.join(self.task_directory,'gulp.in')))
            print("\t",os.path.exists(os.path.join(self.task_directory,'gulp.in')))
            print("poscar:{}".format(
                os.path.join(self.task_directory,self.structure_filename)))
            print("\t",os.path.exists(os.path.join(self.task_directory,self.structure_filename)))
            raise
        gulp_bin = os.environ.get('GULP_BIN')

        cmd_str = '{} < gulp.in > gulp.out 2>/dev/null'.format(gulp_bin)

        orig_dir = os.getcwd()
        os.chdir(self.task_directory)
        try:
            print(cmd_str)
            subprocess.call(cmd_str,shell=True)
        finally:
            os.chdir(orig_dir)

    def write_gulp_input_file(self,filename='gulp.in',poscar='POSCAR'):
        """ writes the gulp input file

        Args:
            filename (str): location to write the gulp input file.  Default is 'gulp.in'
            poscar (str): location of poscar file for structure to read.  
                Default is 'POSCAR'
        """
        self.structure_filename = poscar
        str_out = "opti conp prop phon eigen\n"
        str_out += self.gulp_positions_to_string(poscar)
        str_out += self.potential.gulp_potential_section_to_string(param_dict)
        str_out += self.gulp_phonon_section_to_string()

        try:
            with open(filename,'w') as f:
                f.write(str_out)
        except FileNotFoundError as e:
            print('PYPOSPACK_ERROR: Cannot find filename, {}'.format(filename))
            if os.path.exists(self.task_directory):
                print("\ttask_directory exists")
            else:
                print("\ttask_directory does not exist, {}".format(self.task_directory))
            raise
            
    def gulp_phonon_section_to_string(self,shrink=[8,8,8],kpoints=[10,10,10]):
        str_out = (\
            "shrink\n"
            "{} {} {}\n"
            "kpoints\n"
            "{} {} {}\n"
            "output freq text freq.gulp 12 \n"
            "output phonon text phonon.gulp\n"
            "output osc test phonon.osc\n").format(\
                    shrink[0],shrink[1],shrink[2],
                    kpoints[0],kpoints[1],kpoints[2])

        return str_out

    def gulp_positions_to_string(self,structure='POSCAR'):
        """ returns the gulp position section to create simulation cell

        Args:
            structure (str): the filename of the POSCAR format file.  Default 
                is 'POSCAR'.  If an object which subclasses the 
                pyposmat.crystal.SimulationCell class is passed, this method
                will use this that class instead.

        Returns
            str:
        """

        sim_cell = None
        if isinstance(structure,crystal.SimulationCell):
            sim_cell = crystal.SimulationCell(poscar)
        else:
            sim_cell = vasp.Poscar()
            try:
                sim_cell.read(structure)
            except FileNotFoundError as e:
                msg = 'PYPOSPACK_ERROR: cannot find the POSCAR file:{}'.format(structure)
                raise

        H = sim_cell.H * sim_cell.a0
        str_out = "vectors\n"
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[0,0],H[0,1],H[0,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[1,0],H[1,1],H[1,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[2,0],H[2,1],H[2,2])
        str_out += "fractional\n"
        for s in sim_cell.symbols:
            for a in sim_cell.atomic_basis:
                if a.symbol == s:
                    try:
                        str_out += "{} core {} {} {}\n".format(\
                                s,a.position[0],a.position[1],a.position[2])
                    except:
                        print(s)
                        print(a.symbol)
                        print(a.position)
                        raise
        return str_out

if __name__ == '__main__':
    vasp_filename = 'MgO_NaCl_unit.vasp'
    vasp_input_filename = os.path.join(os.getcwd(),'rsrc',vasp_filename)

    param_dict = {}
    param_dict['chrg_Mg'] = +2.0
    param_dict['chrg_O']  = -2.0
    param_dict['MgMg_A']   = 0.0 
    param_dict['MgMg_rho'] = 0.5
    param_dict['MgMg_C']   = 0.0
    param_dict['MgO_A']    = 821.6
    param_dict['MgO_rho']  = 0.3242
    param_dict['MgO_C']    = 0.0
    param_dict['OO_A']     = 2274.00 
    param_dict['OO_rho']   = 0.1490
    param_dict['OO_C']     = 27.88

    #### TEST INITIALIZE ####
    task_name = 'gulp_test'
    task_directory = os.path.join(\
            os.getcwd(),
            task_name)
    gulp_input_filename = 'gulp.in'
    #task = GulpPhononCalculation(task_name,task_directory)

    #### TEST POTENTIAL POTENTIAL ####
    #print('----- test that buckingham potential provides the right format -----')
    #task.potential = potential.Buckingham(['Mg','O'])
    #task.param_dict = copy.deepcopy(param_dict)
    #print(task.potential.gulp_potential_section_to_string(param_dict))

    #### TEST IF WE CAN CAN WRITE THE INPUT FILE ####
    #gulp_input_filename = os.path.join(task.task_directory,'gulp.in')
    #poscar_input_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    #task.write_gulp_input_file(\
    #        filename=gulp_input_filename,
    #        poscar=vasp_input_filename)

    #### TEST IF WE CAN RUN THE BUCKINGHAM POTENTIAL ####
    #task.structure_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    #task.run()


    task = GulpPhononCalculation(task_name,task_directory)
    task.potential = potential.Buckingham(['Mg','O'])
    task.param_dict = copy.deepcopy(param_dict)
    task.structure_file = os.path.join(\
            task.task_directory,
            vasp_input_filename)
    task.write_gulp_input_file(\
            filename=os.path.join(task.task_directory,gulp_input_filename),
            poscar=vasp_input_filename)
    task.run()
