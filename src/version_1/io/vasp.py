# -*- coding: utf-8 -*-

"""Input and output functions and classes for VASP """
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, shutil, pathlib
import re
import copy
#import pyflamestk.base as base
import pypospack.crystal as crystal
import numpy as np

# *****************************************************************************
# ****    ERROR EXCEPTION HANDLING CLASSES                                 ****
# *****************************************************************************

class VaspIncarError(Exception):
    def __init__(self,*args,**kwargs):
        """Error class for reading/writing VASP INCAR IO issues """
        Exception.__init__(self,*args,**kwargs)

class VaspPoscarError(Exception):
    def __init__(self,*args,**kwargs):
        """Error class for reading/writing VASP POSCAR IO issues """
        Exception.__init__(self,*args,**kwargs)

class VaspPotcarError(Exception):
    def __init__(self,*args,**kwargs):
        """ Error class for reading/writing VASP POTCAR IO issues """
        Exception.__init__(self,*args,**kwargs)

# *****************************************************************************
# ****    SOME CONVENIENCE FUNCTIONS                                       ****
# *****************************************************************************

def read_poscar_file(filename):
    poscar = Poscar()
    poscar.read(filename=filename)
    return poscar

def write_poscar_file(self,obj_poscar,filename='POSCAR'):
    obj_poscar.write(filename)

def read_incar_file(filename):
    incar = Incar()
    incar.read(filename=filename)
    return incar

def make_super_cell(obj, scp):
    sc = base.make_super_cell(copy.deepcopy(obj), list(scp))
    return copy.deepcopy(Poscar(sc))

# *****************************************************************************
# ****     CORE CLASSES                                                    ****
# *****************************************************************************

class VaspSimulation(object):
    """ class for managing vasp simulations

    This class encapsulates are variety of object for a vasp simulation

    Args:
        sim_dir(str): the path of this simulation
        xc(str): the exchange correlation functional

    Attribute:
        poscar (pypospack.io.poscar): class for mananging IO of structures for vasp.
        incar (pypospack.io.incar): class for managing IO of DFT configuration for vasp
        kpoints (pypospack.io.kpoints): defines the the grid for BZ zone integration.
        potcar (pypospack.io.potcar): defines potential file

    """
    def __init__(self,sim_dir,xc='GGA'):
        self.sim_dir = None

        self.is_req_met = False
        self.is_job_pending = True
        self.is_job_initalized = False
        self.is_job_submitted = False
        self.is_job_complete = False
        self.is_job_finished = False

        # vasp input/output filehandlers
        self.poscar = Poscar()
        self.incar = Incar()
        self.kpoints = Kpoints()
        self.potcar = Potcar()
        self.symbols = None

        # set simulation directory
        self.__check_restart_info(sim_dir)
        self.__set_simulation_directory(sim_dir)

    def __check_restart_info(sim_dir):
        if os.path.exists(self.sim_dir):
            pass

    def __set_simulation_directory(self,sim_dir):
        pass

    def read_incar(self,filename='INCAR'):
        self.incar.read(filename=filename)

    def write_incar(self,filename='INCAR'):
        self.incar.write(filename=filename)

    def read_poscar(self,filename='POSCAR'):
        self.poscar.read(filename)

    def write_poscar(self,filename='POSCAR'):
        self.poscar.write(filename)

    def write_kpoints(self,filename='KPOINTS'):
        self.kpoints.write(filename)

    def write_potcar(self,filename='POTCAR'):
        self.potcar.write(filename)

    def read_simulation_cell(self,filename='POSCAR'):
        self._poscar.read(filename)

    def __get_symbols_poscar(self):
        self.symbols = list(self.poscar.symbols)

class VaspMinimizeStructure(object):
    """ class for managing vasp structural minimizations """
    def __init__(self,sim_dir,structure, xc = 'GGA'):
        VaspSimulation.__init__(self)

        if instance(structure,str):
            self.__process_poscar_file(structure)

    def __process_poscar_file():
        pass

class Incar(object):
    def __init__(self, filename="INCAR"):
        """ object for dealing with input and output to VASP via INCAR file

        Args:
        filename (str): the filename of the INCAR file, default:'INCAR'
        """
        self._filename = filename

        self._fmt_section = '# {:*^78}\n'
        self._fmt_arg = '{:<30}! {}\n'
        self._cmt_dict = initialize_incar_comments()

        # default initialization of INCAR file
        self.system = 'automatically generated by pyflamestk'
        self.__init_start_info()
        self.__init_density_of_states()
        self.__init_symmetry()
        self.__init_scf()
        self.__init_spin()
        self.__init_ionic_relaxation()
        self.__init_output()

    def __init_start_info(self):
        self.istart=0
        self.icharg=0

    def __init_density_of_states(self):
        self.ismear=0
        self.sigma=0.2

    def __init_symmetry(self):
        self.isym = 2
        self.symprec = 1e-4
    def __init_scf(self):
        self.ediff = 1e-6 # convergece criteria in eV
        self.nelm = 40 # maximum number of SCF steps
        self.encut = 400 # energy cutoff
        self.prec = 'High'  # set avoid anti-aliasing errors
        self.lreal = 'False' # real space projectors are less accurate
        self.algo = 'Normal' # most robust operator

    def __init_spin(self):
        self.ispin=1
        self.magmom=None

    def __init_ionic_relaxation(self):
        self.ibrion = None
        self.isif = None
        self.ediffg = None
        self.nsw = None
        self.potim = None

    def __init_output(self):
        self.lwave = False
        self.lcharg = False
        self.lvtot = False

    @property
    def is_continue_job(self):
        "bool: True if continuation job"
        if self.istart != 0 or self.icharg!=0:
            return True
        else:
            return False

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, fn):
        self._filename = fn

    def write(self,filename=None):
        """write poscar file

        Args:
        filename (str): the filename of the poscar file,
        """

        if filename is not None:
           self.filename = filename
        f = open(self.filename,'w')
        f.write(self.to_string())
        f.close()

    def read(self,fname=None):
        if fname is not None:
            self.filename = fname
        f = open(self.filename)
        for line in f:
            if line.startswith('#'):
                # ignore comments
                pass
            elif line.strip() == '':
                # ignore blank lines
                pass
            else:
                args = [ line.strip().split('!')[0].split('=')[0].strip(),
                         line.strip().split('!')[0].split('=')[1].strip() ]
                if args[0] == 'ISTART':
                    self.istart = int(args[1])
                elif args[0] == 'ICHARG':
                    self.icharg = int(args[1])
                elif args[0] == 'ISPIN':
                    self.ispin = int(args[1])
                elif args[0] == 'MAGMOM':
                    self.magmom = args[1]
                elif args[0] == 'ISMEAR':
                    self.ismear = int(args[1])
                elif args[0] == 'SIGMA':
                    self.sigma = float(args[1])
                elif args[0] == 'ALGO':
                    self.algo = args[1]
                elif args[0] == 'EDIFF':
                    self.ediff = float(args[1])
                elif args[0] == 'ENCUT':
                    self.encut = float(args[1])
                elif args[0] == 'NELM':
                    self.nelm == int(args[1])
                elif args[0] == 'PREC':
                    self.prec = args[1]
                elif args[0] == 'EDIFFG':
                    self.ediffg = float(args[1])
                elif args[0] == 'IBRION':
                    self.ibrion = int(args[1])
                elif args[0] == 'ISIF':
                    self.isif = int(args[1])
                elif args[0] == 'POTIM':
                    self.potim = float(args[1])
                elif args[0] == 'NSW':
                    self.nsw = int(args[1])
                elif args[0] == 'SYSTEM':
                    self.system = args[1]
                elif args[0] == 'LWAVE':
                    self.lwave = args[1]
                elif args[0] == 'LCHARG':
                    self.lcharg = args[1]
                elif args[0] == 'LVTOT':
                    self.lvtot = args[1]
                elif args[0] == 'LREAL':
                    self.lreal = args[1]
                elif args[0] == 'MAGMOM':
                    self.magmom = args[1]
                elif args[0] == 'ISYM':
                    self.isym = int(args[1])
                elif args[0] == 'SYMPREC':
                    self.symprec = float(args[1])
                else:
                    err_msg = "pypospack does not support tag {}".format(args[0])
                    raise VaspIncarError(err_msg)

    def set_no_ionic_relaxation(self):
        self.ibrion = None
        self.isif = None

    def to_string(self):
        str_out = ''
        str_out += self.__system_information_to_string()
        str_out += self.__start_information_to_string()
        str_out += self.__dos_information_to_string()
        str_out += self.__sym_information_to_string()
        str_out += self.__scf_information_to_string()
        str_out += self.__spin_polarization_to_string()
        str_out += self.__ionic_relaxation_to_string()
        str_out += self.__output_configuration_to_string()
        return str_out

    def __system_information_to_string(self):
        str_out = "SYSTEM = {}\n\n".format(self.system)
        return str_out

    def __start_information_to_string(self):
        fmt = "{} = {}"
        str_out = self._fmt_section.format('STARTING INFORMATION')
        str_out += self._fmt_arg.format(fmt.format('ISTART',self.istart),self._cmt_dict['ISTART'][self.istart])
        str_out += self._fmt_arg.format(fmt.format('ICHARG',self.icharg),self._cmt_dict['ICHARG'][self.icharg])
        str_out += "\n"
        return str_out

    def __dos_information_to_string(self):
        fmt = "{} = {}"
        str_out = self._fmt_section.format('DENSITY OF STATES')
        str_out += self._fmt_arg.format(fmt.format('ISMEAR',self.ismear),self._cmt_dict['ISMEAR'][self.ismear])
        str_out += self._fmt_arg.format(fmt.format('SIGMA',self.sigma),self._cmt_dict['SIGMA'])
        str_out += "\n"
        return str_out

    def __sym_information_to_string(self):
        fmt = "{} = {}"
        str_out = self._fmt_section.format('SYMMETRY')
        str_out += self._fmt_arg.format(fmt.format('ISYM',self.isym),self._cmt_dict['ISYM'][self.isym])
        str_out += self._fmt_arg.format(fmt.format('SYMPREC',self.symprec),self._cmt_dict['SYMPREC'])
        str_out += "\n"
        return str_out


    def __scf_information_to_string(self):
        fmt = "{} = {}"
        if self.algo.startswith('N'):
            self.algo = 'Normal'
        elif self.algo.startswith('V'):
            self.algo = 'VeryFast'
        elif self.algo.startswith('F'):
            self.algo = 'Fast'
        else:
            raise VaspIncarError('Unsupported ALGO value:{}'.format(self.algo))

        if self.lreal == False or self.lreal.startswith('F') or self.lreal == '.FALSE.':
            self.lreal = '.FALSE.'
        elif self.lreal == True or self.lreal.startswith('T') or self.lreal == '.TRUE.':
            self.lreal = '.TRUE.'
        else:
            str_type = type(self.lreal)
            str_value = '{}'.format(self.lreal)
            print(str_type,str_value)
            raise VaspIncarError('Unsupported LREAL value:{}'.format(self.lreal))

        str_out = self._fmt_section.format('ELECTRONIC SCF RELAXATION')
        str_out += self._fmt_arg.format(fmt.format('ALGO',self.algo),self._cmt_dict['ALGO'][self.algo])
        str_out += self._fmt_arg.format(fmt.format('PREC',self.prec),self._cmt_dict['PREC'][self.prec])
        str_out += self._fmt_arg.format(fmt.format('LREAL',self.lreal),self._cmt_dict['LREAL'][self.lreal])
        str_out += self._fmt_arg.format(fmt.format('EDIFF',self.ediff),self._cmt_dict['EDIFF'])
        str_out += self._fmt_arg.format(fmt.format('ENCUT',self.encut),self._cmt_dict['ENCUT'])
        str_out += self._fmt_arg.format(fmt.format('NELM',self.nelm),self._cmt_dict['NELM'])
        str_out += '\n'

        return str_out

    def __spin_polarization_to_string(self):
        fmt = "{} = {}"

        str_out = self._fmt_section.format('SPIN POLARIZATION CONFIGURATION')
        str_out += self._fmt_arg.format(fmt.format('ISPIN',self.ispin),self._cmt_dict['ISPIN'][self.ispin])
        if self.ispin == 2:
            str_out += fmt.format('MAGMOM',self.magmom) + '\n'
        str_out += '\n'

        return str_out

    def __ionic_relaxation_to_string(self):
        fmt = "{} = {}"

        if (self.ibrion is None) and (self.isif is None):
            return ''

        # some default configuration for ibrion
        if self.ibrion is None:
            self.ibrion = 2

        # some default configuration for isif
        if self.isif is None:
            self.isif = 3

        # some default configuration for EDIFFG
        if self.ediffg is None:
            self.ediffg = -0.01 # ev/A - typical
        # some default configuration for POTIM
        if self.potim is None:
            if self.ibrion in [1,2,3]:
                self.potim = 0.5
            if self.ibrion in [5,6,7,8]:
                self.potim = 0.015

        # some default configuration for NSW
        if self.nsw is None:
            self.nsw = 40

        str_out = self._fmt_section.format('IONIC RELAXATION CONFIGURATION')
        str_out += self._fmt_arg.format(fmt.format('IBRION',self.ibrion),self._cmt_dict['IBRION'][self.ibrion])
        str_out += self._fmt_arg.format(fmt.format('ISIF',self.isif),self._cmt_dict['ISIF'][self.isif])
        str_out += self._fmt_arg.format(fmt.format('POTIM',self.potim),self._cmt_dict['POTIM'])
        str_out += self._fmt_arg.format(fmt.format('NSW',self.nsw),self._cmt_dict['NSW'])

        if self.ediffg < 0:
            str_out += self._fmt_arg.format(fmt.format('EDIFFG',self.ediffg),'force convergence requirements in ev A')
        else:
            str_out += self._fmt_arg.format(fmt.format('EDIFGG',self.ediffg),'energy convergence in eV')
        str_out += '\n'

        return str_out

    def __output_configuration_to_string(self):
        fmt = "{} = {}"

        if self.lwave in [True,'T','.TRUE.']:
            self.lwave = '.TRUE.'
        elif self.lwave in [False,'F','.FALSE.']:
            self.lwave = '.FALSE.'
        else:
            msg = "LWAVE tag cannot be {}({})".format(type(self.lwave),self.lwave)
            raise VaspIncarException(msg)

        if self.lcharg in [True,'T','.TRUE.']:
            self.lcharg = '.TRUE.'
        elif self.lcharg in [False,'F','.FALSE.']:
            self.lcharg = '.FALSE.'
        else:
            msg = "LCHARG tag cannot be {}({})".format(type(self.lcharg),self.lcharg)
            raise VaspIncarException(msg)

        if self.lvtot is [True,'T','.TRUE.']:
            self.lvtot = '.TRUE.'
        elif self.lvtot in [False,'T','.FALSE.']:
            self.lvtot= '.FALSE.'
        else:
            msg = "LVTOT tag cannot be {}({})".format(type(self.lwave),self.lvtot)
            raise VaspIncarException(msg)

        str_out = self._fmt_section.format('OUTPUT CONFIGURATION')
        str_out += self._fmt_arg.format(fmt.format('LWAVE',self.lwave),self._cmt_dict['LWAVE'][self.lwave])
        str_out += self._fmt_arg.format(fmt.format('LCHARG',self.lcharg),self._cmt_dict['LCHARG'][self.lcharg])
        str_out += self._fmt_arg.format(fmt.format('LVTOT',self.lvtot),self._cmt_dict['LVTOT'][self.lvtot])
        str_out += "\n"
        return str_out

class Outcar(object):
    """ object to manage IO to VASP through the OUTCAR file

    Many of the attributes in this class are initialized to None at
    instanciation.  They are processed when the OUTCAR file is read.

    Args:
        filename (str): filename of the OUTCAR file

    Attributes:
        total_energy (float): total energy of the structure in eV.
        elastic_tensor (nd.array): the components of the elastic tensor.  The
            components of the array start at :math:`0`.  So the :math:`c_{11}`
            component would be accessed by elastic_tensor[0,0]
        phonon_eig_val (nd.array): the eigenvalues of the phonons determined
            by lattice dynamics.
        phonon_eig_vec (nd.array): the eigenvectors of the phonons determined
            by lattice dynamics.
    """
    def __init__(self, filename="OUTCAR"):
        self.filename = filename
        self.total_energy = None
        self.encut = None
        self.entropy = None
        self.elastic_tensor = None
        self.phonon_eig_val = None # phonon eigenvalues
        self.phonon_eig_vec = None # phonon eigenvectors

    def read(self,filename=None):
        """ read a VASP outcar file and parse for information

        Args:
            filename (str): filename
            encut (float): energy cutoff for this simulation
            total_energy (float): total energy for this simulation
        """

        if filename is not None:
            self.filename = filename
        with open(self.filename) as f:
            for line in f:
                # check if free energy line
                if "TOTEN" in line:
                    try:
                        E = line.strip().split('=')[1].strip().split(' ')[0]
                        E = float(E)
                        self.total_energy = E
                    except ValueError as e:
                        if type(self.total_energy) is not float:
                            pass
                    except IndexError as e:
                        if type(self.total_energy) is not float:
                            pass
                elif "ENCUT" in line:
                    E = line.strip().split('=')[1].strip().split(' ')[0]
                    E = float(E)
                    self.encut = E
                elif "EENTRO" in line:
                    try:
                        E = line.strip().split('=')[1].strip()
                        E = float(E)
                        self.entropy = E
                    except ValueError as e:
                        if type(self.entropy) is not float:
                            pass
                    except IndexError as e:
                        if type(self.entropy) is not float:
                            pass
                else:
                    pass

    def __string__():
        print("total_energy[eV] = {}".format(self._ecoh_per_structure))

    def get_phonons(self):
        pass

    def get_time_of_calculation(self):
        pass

    def get_number_of_atoms(self):
        pass

    def get_energy(self, fname="OUTCAR"):
        self.ecoh_per_structure = None
        self.ecoh_per_atom = None

class Potcar(object):
    """ object to deal with POTCAR files for VASP simulations

    Args:
        symbols (:obj:`list` of :obj:`str`): a list of chemical symbols
        filename (str): the filename of the POTCAR
        xc (str): the exchange correlation model.  Default is 'GGA'

    Attributes:
        potcar_dir (str): the path of the VASP potential directory.
        potcars (dict): a dictionary of key-value pairs of symbols to potcar
            type where the key value is the ISO symbol, and the value is the
            the psuedopotnetial type.

    If an :obj:`str` is passed into symbols, then the string will be
    cast into a :obj:`list` of :obj:`str`

    A directory of potentials is not included in pypospack because it is
    protected by copy right.  If the potcar_dir is not set, then pypospack
    will use the environment variables:

    >>> export VASP_LDA_DIR=$(cd ~/opt/vasp/potentials/LDA;pwd)
    >>> export VASP_GGA_DIR=$(cd ~/opt/vasp/potentials/GGA;pwd)

    """
    def __init__(self, symbols = None,
                       filename= None,
                       xc = 'GGA',
                       potcars = None):
        self.potcar_dir = None

        self.filename = None
        self.symbols = None
        self.potcars = None

        if filename is not None:
            self.filename = filename

        if isinstance(symbols,str):
            self.symbols = [symbols]
        else:
            try:
                self.symbols = [s for s in symbols]
            except:
                msg = 'symbols type{}={}'.format(type(symbols),symbols)

        if potcars is not None:
            self.potcars = potcars.copy()
        self.xc = xc
        self.encut_min = []
        self.encut_max = []
        self.models    = []
        self.exch      = []

    def read(self, fname = "POTCAR"):
        # initialize arrays
        self.symbols   = []
        self.encut_min = []
        self.enmin_max = []
        self.models = []
        self.xc      = []

        with open(fname,'r') as f:
            for line in f:
                line = line.strip()
                if 'TITEL' in line:
                    symbol = line.split('=')[1].strip().split(' ')[1]
                    self.symbols.append(symbol)
                elif 'ENMIN' in line:
                    obj_re = re.match('ENMAX  =  (.*); ENMIN  =  (.*) eV.*',line,re.M|re.I)
                    enmax = float(obj_re.group(1))
                    enmin = float(obj_re.group(2))
                    self.encut_min.append(enmin)
                    self.encut_max.append(enmax)
                elif 'LEXCH' in line:
                    xc = line.split('=')[1].strip()
                    self.xc.append(xc)
                elif "VRHFIN" in line:
                    m = line.split('=')[1].strip()
                    self.models.append(m)

    def write(self,filename = None, src = None, potcars = None):
        """ write POTCAR file

        Raises:
            pypospack.io.vasp.VaspPotcarError
        """
        if filename is not None:
            self.filename = filename
        else:
            if self.filename is None:
                raise TypeError('cannot find suitable filename')

        #  if a source potcar is given just copy the file
        if src is not None:
            __write_potcar_by_copy(src_fn=src,dst_fn=self.filename)

        self.__set_xc_directory()

        #  get full path of potcars
        if potcars is not None:
            self.potcars = potcars.copy()
        else:
            self.__find_potcar_files()


        with open(self.filename,'w') as f_out:
            for s in self.symbols:
                with open(self.potcars[s],'r') as f_in:
                    # avoid reading large files into memory
                    chunk_sz = 1024*1024*10 # 10MB
                    shutil.copyfileobj(f_in,f_out, chunk_sz)

    def __write_potcar_by_copy(self,src_fn,dst_fn):
        with open(dst_fn) as f_out:
            with open(src_fn) as f_in:
                chunk_sz = 1024*1024*10
                shutil.copyfileobj(f_in,f_out,chunk_sz)

    def __find_potcar_files(self):
        # a more intelligent way of selecting potcars based upon VASP
        # recommendations should be implemented
        # http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html
        self.potcars = {} # initialize
        for s in self.symbols:
            # try VASP_XC_DIR/symbol/POTCAR
            if pathlib.Path(os.path.join(self.potcar_dir,s,'POTCAR')).is_file():
                self.potcars[s] = os.path.join(self.potcar_dir,s,'POTCAR')
            # try VASP_XC_DIR/symbol_new/POTCAR
            elif pathlib.Path(os.path.join(self.potcar_dir,"{}_new".format(s),'POTCAR')).is_file():
                self.potcars[s] = os.path.join(self.potcar_dir,"{}_new".format(s),'POTCAR')
            # try VASP_XC_DIR/symbols_h/POTCAR
            elif pathlib.Path(os.path.join(self.potcar_dir,"{}_h".format(s),'POTCAR')).is_file():
                self.potcars[s] = os.path.join(self.potcar_dir,"{}_h".format(s),'POTCAR')
            else:
                print(os.path.join(self.potcar_dir,"{}_new".format(s),'POTCAR'))
                msg = 'cannot find a POTCAR file for {}.{}\n'.format(self.xc,s)
                msg += 'potcar_dir:{}'.format(self.potcar_dir)
                raise VaspPotcarError(msg)

    def __set_xc_directory(self):
        """ sets the directory of the exchange correlation functional

        the exchange correlation functional is set to either 'LDA' or 'GGA'
        based upon environment variables

        Raises:
            (pypospack.io.vasp.VaspPotcarError): if the directory is not set
                 as an environment variable.
        """
        if self.xc == 'GGA':
            try:
                self.potcar_dir = os.environ['VASP_GGA_DIR']
            except KeyError:
                msg = ('need to set environment variable VASP_GGA_DIR to the '
                       'location of the VASP GGA-PBE potential files' )
                raise VaspPotcarError(msg)
        elif self.xc == 'LDA':
            try:
                self.potcar_dir = os.environ['VASP_LDA_DIR']
            except KeyError:
                msg = ('need to set environment variable VASP_LDA_DIR to the '
                       'location of the VASP LDA-CA potential files' )

    def __str__(self):
        header_row   = "symbol enmin enmax xc\n"
        format_row   = "{}({}) {:10.6f} {:10.6f} {}\n"

        n_atoms = len(self._symbols)
        str_out      = header_row
        for i in range(n_atoms):
            str_out += format_row.format(self._symbols[i],
                                         self._models[i],
                                         self._encut_min[i],
                                         self._encut_max[i],
                                         self._xc[i])
        return str_out

class Kpoints(object):
    """ KPOINTS file

    Args:
        filename(str): the filename to either read from or to write to.  Default is KPOINTS
        comment(str): a comment to use for the KPOINTS file.  Default is 'Automatic Mesh'
        n_kpoints(int): number of kpoints for automatic mesh.  Default is 0.
        mesh_type(str): mesh type to create.  Default is 'Monkhorst-Pack'
        mesh_size(:obj:`list` of :obj:`int`): defines the density of the sampling in the
            Brillouin zone of the reciprocal lattice.
        mesh_shift(:obj:`list` of :obj:`int`): defines the shift of the k-kpoint mesh

    Attributes:
        filename(str): the filename to either read from or to write to.  Default is KPOINTS
        comment(str): a comment to use for the KPOINTS file.  Default is 'Automatic Mesh'
        n_kpoints(int): number of kpoints for automatic mesh.  Default is 0.
        mesh_type(str): mesh type to create.  Default is 'Monkhorst-Pack'
        mesh_size(:obj:`list` of :obj:`int`): defines the density of the sampling in the
            Brillouin zone of the reciprocal lattice.
        mesh_shift(:obj:`list` of :obj:`int`): defines the shift of the k-kpoint mesh
    """
    def __init__(self,
                 filename="KPOINTS",
                 comment="Automatic Mesh",
                 n_kpoints = 0,
                 mesh_type = "Monkhorst-Pack",
                 mesh_size = [4,4,4],
                 mesh_shift = [0,0,0]):

        self.filename = filename
        self.comment = comment
        self.n_kpoints = n_kpoints
        self.mesh_type = mesh_type
        self.mesh_size = [i for i in mesh_size]
        self.mesh_shift = [i for i in mesh_shift]


    def to_string(self):
        str_out = "{}\n".format(self.comment)
        str_out += "{}\n".format(self.n_kpoints)
        str_out += "{}\n".format(self.mesh_type)

        str_out += "{:3d} {:3d} {:3d}\n".format(self.mesh_size[0],
                                                self.mesh_size[1],
                                                self.mesh_size[2])
        str_out += "{:3d} {:3d} {:3d}\n".format(self.mesh_shift[0],
                                               self.mesh_shift[1],
                                               self.mesh_shift[2])
        return str_out

    def write(self,filename = None):
        if filename is not None:
            self.filename = filename
        f = open(self.filename, 'w')
        f.write(self.to_string())
        f.close()

    def read(self,filename = None):
        if filename is not None:
            self.filename = filename

class Poscar(crystal.SimulationCell):
    """ POSCAR structure file for VASP

        Args:
            obj_cell (pypospack.crystal.SimulationCell):

        Attributes:
            filename (str): filename (either to be written to or to write from)
            comment (str): comment used in the first line of the POSCAR file
    """

    def __init__(self,obj_cell=None):
        crystal.SimulationCell.__init__(self,obj_cell)
        self.filename = 'POSCAR'
        self.comment = 'Automatically generated by pypospack'

    def read(self, filename=None):
        """ read the POSCAR file

        Args:
            filename (str): the file name of the POSCAR file
        """
        if filename is not None:
            self.filename = filename

        try:
            f = open(self.filename, 'r')
        except FileNotFoundError as e:
            str_out = "\n".join(
                [
                "cwd={}".format(os.getcwd()),
                "filename={}".format(self.filename)
                ])
            raise
        # read structure comment
        self.comment = f.readline().strip()

        # read lattice parameter
        try:
            line = f.readline()
            self.a0 = float(line)
        except ValueError as e:
            msg_err = "Cannot read the lattice parameter from the POSCAR file\n"
            msg_err += "filename:{}\n".format(self.filename)
            msg_err += "line({}):\'{}\'".format(
                str(type(line)),
                line)
            print(msg_err)
            raise ValueError(line)

        h_matrix = np.zeros(shape=[3,3])
        for i in range(3):
            h_row = f.readline().strip().split()
            h_row = np.array([float(val) for val in h_row])
            h_matrix[i,:] = h_row
        self.H = h_matrix.copy()

        # read symbols
        symbols = f.readline().strip().split()

        # number of atoms per symbol
        n_atoms_per_symbol = {}
        line = f.readline().strip().split()
        for i,n in enumerate(line):
            s = symbols[i]
            n_atoms_per_symbol[s] = int(n)

        # read in coordinate type
        line = f.readline()
        coordinate_style = line[0].upper()
        if coordinate_style == "D":
            self.coordinate_style = "Direct"
        elif coordinate_style == "C":
            self.coordinate_style = "Cartesian"
        else:
            msg = "unable to determine the coordinate style {}".format(line)
            raise VaspPoscarError(msg)

        # read in atomic positions
        for s in symbols:
            n_atoms = n_atoms_per_symbol[s]
            for i_atom in range(n_atoms):
                line = f.readline().strip().split()
                position = [float(line[i]) for i in range(3)]
                try:
                    self.add_atom(s,position)
                except:
                    raise
    def write(self, filename=None):
        """ write poscar file """
        if filename is None:
            filename = self.filename
        self.filename = filename

        str_poscar  = self.comment + "\n"
        str_poscar += "{:10.6}\n".format(self.a0)

        # string for h-matrix
        for i in range(3):
            h_row_template = "{:10.6f} {:10.6f} {:10.6f}\n"
            str_poscar += h_row_template.format(self.H[i,0],
                                                self.H[i,1],
                                                self.H[i,2])

        sym_list = self.symbols
        str_atomlist = ""
        str_atomnum  = ""
        for sym in sym_list:
            nAtoms = self.get_number_of_atoms(sym)
            str_atomlist   += " " + sym
            str_atomnum    += " " + str(nAtoms)
        str_atomlist   += "\n"
        str_atomnum    += "\n"

        str_poscar += str_atomlist
        str_poscar += str_atomnum
        str_poscar += "Direct\n"

        for symbol in self.symbols:
            for i, atom in enumerate(self.atomic_basis):
                if symbol == atom.symbol:
                    pos_template = "{:10.6f} {:10.6f} {:10.6f}\n"
                    str_position = pos_template.format(atom.position[0],
                                                       atom.position[1],
                                                       atom.position[2])
                    str_poscar += str_position

        f = open(self.filename, 'w')
        f.write(str_poscar)
        f.close()

# *****************************************************************************
# ****    SOME HELPER FUNCTIONS
# *****************************************************************************

def get_incar_istart_comments(option):
    options = {
        0:'begin from scratch',
        1:'continuation job, constant energy cutoff',
        2"'continuation job, constant basis set'
    }
    return options[option]

def get_incar_isym_comments(option):
    options = {
        -1:'symmetry off',
        0:'symmetry on',
        1:'symmetry on',
        2:'symmetry on, efficient symmetrization',
        3:'symmetry on, only forces and stress tensor'
    }
    return options[option]

def get_incar_icharg_comments(option):
    options = {
        0:'Calculate charge density from initial wave functions.',
        1:'Read the charge density from file CHGCAR',
        2:'Take superposition of atomic charge densities'
    }
    return options[option]

def get_incar_ismear_comments(option):
    options = {
        -5:'tetrahedron method with Blochl corrections',
        -4:'tetrahedron method',
        0:'method of Gaussian smearing',
        1:'method of Methfessel-Paxton order 1',
        2:'method of Methfessel-Paxton order 2'
    }
    return  options[option]

def get_incar_sigma_comment():
    return 'width of the smearing in eV.'

def initialize_incar_comments():
    comments_dict = {}
    comments_dict["SYSTEM"] = ''
    comments_dict["ISTART"] = {}
    comments_dict["ISTART"][0] = 'begin from scratch'
    comments_dict["ISTART"][1] = 'continuation job, constant energy cutoff'
    comments_dict["ISTART"][2] = 'continuation job, constant basis set'

    comments_dict["ISYM"] = {}
    comments_dict["ISYM"][-1] = 'symmetry off'
    comments_dict["ISYM"][0] = 'symmetry on'
    comments_dict["ISYM"][1] = 'symmetry on'
    comments_dict["ISYM"][2] = 'symmetry on, efficient symmetrization'
    comments_dict["ISYM"][3] = 'symmetry on, only forces and stress tensor'

    comments_dict["SYMPREC"] = 'determines how accurate positions must be'

    comments_dict["ICHARG"] = {}
    comments_dict["ICHARG"][0] = 'Calculate charge density from initial wave functions.'
    comments_dict["ICHARG"][1] = 'Read the charge density from file CHGCAR'
    comments_dict["ICHARG"][2] = 'Take superposition of atomic charge densities'

    comments_dict["ISMEAR"] = {}
    comments_dict["ISMEAR"][-5] = 'tetrahedron method with Blochl corrections'
    comments_dict["ISMEAR"][-4] = 'tetrahedron method'
    comments_dict["ISMEAR"][0] = 'method of Gaussian smearing'
    comments_dict["ISMEAR"][1] = 'method of Methfessel-Paxton order 1'
    comments_dict["ISMEAR"][2] = 'method of Methfessel-Paxton order 2'

    comments_dict["SIGMA"] = 'width of the smearing in eV.'

    comments_dict["NELM"] = 'maximum number of electronic SC'
    comments_dict["ENCUT"] = 'Cut-off energy for plane wave basis set in eV'
    comments_dict["EDIFF"] = 'convergence condition for SC-loop in eV'

    comments_dict["PREC"] = {}
    comments_dict["PREC"]['Accurate'] = 'avoid wrap around errors'
    comments_dict["PREC"]['High'] = 'avoid wrap around errors'

    comments_dict["ALGO"] = {}
    comments_dict["ALGO"]["Normal"] = 'blocked Davidson iteration scheme'
    comments_dict["ALGO"]["VeryFast"] = 'RMM-DIIS'
    comments_dict["ALGO"]["Fast"] = 'blocked Davidson, followed by RMM_DIIS'

    comments_dict["LREAL"] = {}
    comments_dict["LREAL"]['.FALSE.'] = 'projection done in reciprocal space'
    comments_dict["LREAL"]['On'] = 'method of King-Smith, et al. Phys. Rev B 44, 13063 (1991).'
    comments_dict["LREAL"]['Auto'] = 'unpublished method of G. Kresse'

    comments_dict["ISPIN"] = {}
    comments_dict["ISPIN"][1] = 'non-spin polarized calculations'
    comments_dict["ISPIN"][2] = 'spin polarized calculations'

    comments_dict["IBRION"] = {}
    comments_dict["IBRION"][0] = 'molecular dynamics'
    comments_dict["IBRION"][1] = 'ionic relaxation by RMM-DIIS'
    comments_dict["IBRION"][2] = 'ionic relaxation by CG'
    comments_dict["IBRION"][3] = 'ionic relaxation by damped MD'
    comments_dict["IBRION"][5] = 'phonons, by frozen ion, without symmetry'
    comments_dict["IBRION"][6] = 'phonons, by frozen ion, with symmetry'
    comments_dict["IBRION"][7] = 'phonons, by perturbation theory, no symmetry'
    comments_dict["IBRION"][8] = 'phonons, by perturbtion theory, with symmetry'

    comments_dict["ISIF"] = {}
    comments_dict["ISIF"][2] = 'relaxation, ions=T, cellshape=F, cellvolume=F'
    comments_dict["ISIF"][3] = 'relaxation, ions=T, cellshape=T, cellvolume=T'
    comments_dict["ISIF"][4] = 'relaxation, ions=T, cellshape=T, cellvolume=F'
    comments_dict["ISIF"][5] = 'relaxation, ions=F, cellshape=T, cellvolume=F'
    comments_dict["ISIF"][6] = 'relaxation, ions=F, cellshape=T, cellvolume=T'
    comments_dict["ISIF"][7] = 'relaxation, ions=F, cellshape=F, cellvolume=T'

    comments_dict["POTIM"] = 'scaling factor in relaxation'
    comments_dict["NSW"] = 'maximum number of ionic relaxation steps'
    comments_dict["LWAVE"] = {}
    comments_dict["LWAVE"]['.TRUE.'] = 'write WAVECAR'
    comments_dict["LWAVE"]['.FALSE.'] = 'do not write WAVECAR'

    comments_dict["LCHARG"] = {}
    comments_dict["LCHARG"]['.TRUE.'] = 'write CHGCAR, write CHG'
    comments_dict["LCHARG"]['.FALSE.'] = 'no CHGCAR, no CHG'

    comments_dict['LVTOT'] = {}
    comments_dict['LVTOT']['.TRUE.'] = 'write LOCPOT'
    comments_dict['LVTOT']['.FALSE.'] = 'no LOCPOT'
    return comments_dict

# standard recommendations
potcar_std = { 'H':'H',
               'He':'He',
               'Li':'Li_sv',
               'Be':'Be',
               'B':'B',
               'C':'C'}
