import os,copy,shutil,subprocess
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

from pypospack.task.gulp import GulpSimulationError, GulpSimulation

class GulpPhononCalculation(GulpSimulation):
    def __init__(
            self,
            task_name,
            task_directory,
            structure_filename,
            restart=False
        ):
        # initialize the parent class
        GulpSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                structure_filename=structure_filename,
                )

        self.shrink=[8,8,8]
        self.kpoints=[10,10,10]

    def setup_task_directory(self):
        self.gulp_input_filename = os.path.join(
                self.task_directory,'gulp.in')
        self.write_gulp_input_file(
                filename=self.gulp_input_filename,
                poscar=self.structure_filename)

    def run(self):
        cmd_str = '{} < gulp.in > gulp.out 2>/dev/null'.format(self.gulp_bin)

        original_directory = os.getcwd()
        os.chdir(self.task_directory)

        try:
            subprocess(cmd_str,shell=True)
        finally:
            os.chdir(original_directory)

    def write_gulp_input_file(self,filename=None,structure_filename=None):
        """
        Args:
            filename (str): location to write the gulp input file.
            poscar (str): location of poscar file for structure to read.  
        """

        if filename is not None:
            self.gulp_input_filename = filename
        if structure_filename is not None:
            self.structure_filename=structure_filename

        str_out = "opti conp prop phon eigen\n"
        str_out += self.get_structure_section(filename)
        str_out += self.get_potential_section(parameters)
        str_out += self.get_phonon_section()

        with open(self.gulp_input_filename,'w') as f:
            f.write(str_out)

    def get_structure_section(self,filename=None):

        if filename is not None:
            self.structure_filename = filename

        # read in the poscar file
        self.structure = vasp.Poscar()
        self.structure.read(self.structure_filename)

        H = sim_cell.H * sim_cell.a0
        str_out = "vectors\n"
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[0,0],H[0,1],H[0,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[1,0],H[1,1],H[1,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[2,0],H[2,1],H[2,2])
        str_out += "fractional\n"
    

    def get_phonon_section(self):

        str_phonon = (
                "shrink\n"
                "{shrink1} {shrink2} {shrink3}\n"
                "kpoints\n"
                "{k1} {k2} {k3}\n"
                "output freq text freq.gulp 12 \n"
                "output phonon text phonon.gulp\n"
                "output osc test phonon.osc\n"
                ).format(
                        shrink1=self.shrink[0],
                        shrink2=self.shrink[1],
                        shirnk3=self.shrink[2],
                        k1=self.kpoints[0],
                        k2=self.kpoints[1],
                        k3=self.points[2])
        return str_out
