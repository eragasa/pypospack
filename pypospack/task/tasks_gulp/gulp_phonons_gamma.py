import os,copy,shutil,subprocess,yaml
from collections import OrderedDict
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

from pypospack.task.gulp import GulpSimulationError, GulpSimulation

class GulpGammaPointPhonons(GulpSimulation):
    """ GULP gamma point phonon calculations.

    This class sets up, runs and processes the phonon calculation at the gamma point only 
    using gulp.  The simulation cell is first structurally minimized before doing the phonon
    calculation.

    Args:
        task_name (str): the name of the task
        task_directory (str): the directory of the task
        debug(bool): by default set to false, if set to True outputs debug
            information to standard out
    """

    def __init__(
            self,
            task_name,
            task_directory,
            structure_filename='POSCAR',
            restart=False,
            fullauto=False,
            debug=False
        ):

        _task_type = 'gulp_gamma_phonons'
        # initialize the base class
        GulpSimulation.__init__(
                self,
                task_name=task_name,
                task_directory=task_directory,
                task_type=_task_type,
                structure_filename=structure_filename,
                restart=restart)

        # set additional attributes
        self.is_debug = debug

        self.n_phonons = None

    
    def postprocess(self,configuration=None):
        GulpSimulation.postprocess(self)

    def on_init(self,configuration=None):
        GulpSimulation.on_init(self,configuration=configuration)

    def on_config(self,configuration=None):
        GulpSimulation.on_config(self,configuration=configuration)

    def on_ready(self, configuration=None, results=None):
        GulpSimulation.on_ready(self,
                configuration=configuration,
                results=results)

    def on_post(self,configuration=None):
        self.__get_phonon_frequencies_from_phonon_frequency_file()
        GulpSimulation.on_post(self,configuration=configuration)

    def on_running(self,configuration=None):
        GulpSimulation.on_running(self,configuration=None)

    def on_finished(self, configuration=None):
        pass

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
    
    def get_gulpinputfile_phononsection(self):

        str_phonon = (
                "output freq text freq.gulp\n"
                "output phonon text phonon.gulp\n"
                "output osc test phonon.osc\n"
                )
        return str_phonon

    def __get_phonon_frequencies_from_phonon_frequency_file(self,phonon_freq_fn=None): 

        # set filename
        if phonon_freq_fn is None:
            _fn = os.path.join(
                self.task_directory,
                'freq.gulp'
                )
        else:
            _fn = phonon_freq_fn

        self.results = OrderedDict()

        #read the filename
        with open(_fn,'r') as _f:
            lines = _f.readlines()
        
        #read the filename
        freqs = [float(line) for line in lines]

        _task_name = self.task_name
        for idx,freq in enumerate(freqs):
            _result_name = 'gamma_phonon_{}'.format(idx+1)

            _rk = '{}.{}'.format(_task_name,_result_name)
            _rv = freq
            
            self.results[_rk] = _rv
        
        self.n_phonons = len(freqs)

