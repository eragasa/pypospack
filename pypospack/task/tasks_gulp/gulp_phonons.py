import os,copy,shutil,subprocess,yaml
from collections import OrderedDict
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
        self.shrink=[8,8,8]
        self.kpoints=[10,10,10]
        # initialize the parent class
        GulpSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                structure_filename=structure_filename,
                )

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
        str_out += self.get_gulpinputfile_structuresection()
        str_out += self.get_gulpinputfile_potentialsection()
        str_out += self.get_gulpinputfile_phononsection()

        gulp_input_filename = os.path.join(
                self.task_directory,
                self.gulp_input_filename)
        with open(gulp_input_filename,'w') as f:
            f.write(str_out)

    def get_conditions_post(self):
        GulpSimulation.get_conditions_post(self)
        self.conditions_POST['is_freq_file_exists'] \
                = os.path.isfile(os.path.join(
                    self.task_directory,
                    'freq.gulp'))
        self.conditions_POST['is_dens_file_exists'] \
                = os.path.isfile(os.path.join(
                    self.task_directory,
                    'phonon.gulp.dens'))

    def get_gulpinputfile_phononsection(self):

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
                        shrink3=self.shrink[2],
                        k1=self.kpoints[0],
                        k2=self.kpoints[1],
                        k3=self.kpoints[2])
        return str_phonon

class GulpGammaPointPhonons(GulpPhononCalculation):

    def get_gulpinputfile_phononsection(self):

        str_phonon = (
                "output freq text freq.gulp 12 \n"
                "output phonon text phonon.gulp\n"
                "output osc test phonon.osc\n"
                )
        return str_phonon

    def on_post(self):
        self.results = OrderedDict()
        freq_filename = os.path.join(
                self.task_directory,
                'freq.gulp')

        lines = None
        with open(freq_filename,'r') as f:
            lines = f.readlines()
        
        freqs = [float(line) for line in lines]
        
        for idx,freq in enumerate(freqs):
            key = "{}.freq.{}".format(self.task_name,idx+1)
            self.results[key] = freq
        
        _results_filename = self.results_filename
        
        with open(_results_filename,'w') as f:
            yaml.dump(self.results,f,default_flow_style=True)
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_finished(self):
        self.cleanup()

    def cleanup(self):
        print(self.results_filename)
        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
