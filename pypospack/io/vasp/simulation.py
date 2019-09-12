import os
from pypospack.io.vasp import Poscar
from pypospack.io.vasp import Incar
from pypospack.io.vasp import Kpoints
from pypospack.io.vasp import Potcar
from pypospack.io.slurm import SlurmSubmissionScript

class VaspSimulation():

    def __init__(self, simulation_path="."):
        self.simulation_path = simulation_path
        self.poscar = Poscar()
        self.incar = Incar()
        self.kpoints = Kpoints()
        self.potcar = Potcar()
        self.submission_script = SlurmSubmissionScript()
    @property
    def symbols(self):
        return self.incar.symbols

    @property
    def xc(self):
        return self.potcar.xc

    @xc.setter
    def xc(self, xc):
        self.potcar.xc = xc

    def read_incar(self,filename=None):
        if filename is None:
            filename_ = os.path.join(self.simulation_path, "INCAR")
        else:
            filename_ = filename

        self.incar.read(filename=filename_)

    def write_incar(self,filename=None):
        if filename is None:
            filename_ = os.path.join(self.simulation_path, "INCAR")
        else:
            filename_ = filename

        self.incar.write(filename=filename_)

    def write_poscar(self, filename=None):
        if filename is None:
            filename_ = os.path.join(self.simulation_path, "POSCAR")
        else:
            filename_ = filename

        self.poscar.write(filename=filename_)

    def write_potcar(self, xc='GGA', filename=None):
        self.xc = xc

        if filename is None:
            filename_ = os.path.join(self.simulation_path, 'POTCAR')
        else:
            filename_ = filename

        self.potcar.symbols = self.symbols
        self.potcar.write(filename=filename_)

    def write_kpoints(self, filename=None):
        if filename is None:
            filename_ = os.path.join(self.simulation_path,"KPOINTS")
        else:
            fiename_ = filename

        self.kpoints.write(filename_)

    def write_submission_script(self, filename=None):
        if filename is None:
            filename_ = os.path.join(self.simulation_path,"runjob.slurm")
        else:
            filename_ = filename

        self.simulation_script.write(filename=filename_)

    def write(self, path, is_clobber=False):
        self.simulation_path = path

        try:
            os.path.mkdir(path)
        except FileExistsError as e:
            if is_clobber:
                shutil.rmthree(path)
                os.path.mkdir(path)
            else:
                raise

        self.write_poscar()
        self.write_potcar()
        self.write_incar()
        self.write_kpoints()
        self.write_submission_script()
