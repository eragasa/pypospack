class Poscar(object):
    pass

class Incar(object):
    algo_types = ['Normal','Fast','Very Fast']

    def __init__(self,filename='INCAR'):
        self.str_format = "{:40}! {}"

        self.istart = 0
        self.encut = None
        self.nelm = None
        self.ediff = 1e-6
        self.ismear = None
        self.ismear_order = None
        self.sigma = None
        self.potim = None
        self.isif = None
        self.algo = 'Normal'
        self.ediffg = None

        self._is_neb = False
        self._is_ionic_relaxation = False

    def write(self,filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError

    def set_ionic_relaxation(self):
        self._is_ionic_relaxation = True
        raise NotImplementedError

    def set_nudged_elastic_band(self,\
            fname_neb_init,
            fname_neb_final,
            n_neb_images,):
        self._is_neb = True

        self.neb_filename_image_initial
        self.neb_filename_image_final
        self.neb_n_images
        raise NotImplementedError

    def set_density_of_states(self):
        raise NotImplementedError

    def _get_density_of_states_string(self):
        raise NotImplementedError
    def read(self,filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError
kpoint_type = ['MonkhorstPack','Gamma']

class Kpoints(object);
    def __init__(self,filename='KPOINTS'):
        self._filename = filename
        self._generation_type = Auto
        self._generation_flag = 0
        self._kpoints = [8,8,8]
        self._shift = [0,0,0]

    def read(self,filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError
    def write(self,filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError

class Potcar(object):
    def __init__(self,filename='POTCAR'):
        self._filename = filename

    def read(self, filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError

    def write(self, filename=None):
        if filename is not None:
            self._filename = filename
        raise NotImplementedError

class Outcar(object):
    def __init__(self,filename='OUTCAR'):
        self._filename = filename
