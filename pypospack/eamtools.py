import copy
import numpy as np
import scipy.constants as spc
from collections import OrderedDict
import pypospack.io.filesystem as filesystem
# import pyflamestk.potentials
#import pyflamestk.tableofelements as toe

def get_eam_radii_vector(r_max,N_r):
    assert type(r_max) == float
    assert type(N_r) == float

    if N_r%5 != 0:
        raise ValueError('N_r must be divisible by 5')

    _r = r_max * np.linspace(1,N_r,N_r)/N_r
    _dr = r_max / N_r

    return _r, _dr

class EamCurveFitter(object):
    
    def __init__(self, element_names):

        self.elements = {}
        self.pair_potentials = {}
        for element_name in element_names:
          print(element_name)
          self.elements[element_name] = Element("Ni")

        # initialize varaibles
        self.embedding_energy = {}
        self.pair_potential = {}
        self.electron_density = {}
        
        print("initializing variables...")

        for el_i in self.elements.keys():
          print('electron_density: rho_{}'.format(el_i))
          print('embedding_energy: F_{}'.format(el_i))
          self.embedding_energy[el_i] = []
          self.electron_density[el_i] = []

        for el_i in self.elements.keys():
          for el_j in self.elements.keys():
            print('pair_potential: p_{}{}'.format(el_i,el_j))
            self.pair_potentials['{}{}'.format(el_i,el_j)] = []

class EamFuncflFile(object):
    pass

class EamSetflFile(object):
    """

    This class wites a setfl file.
   
    There are three components to the EAM model: the embedding energy, the 
    electron density, and interatomic pair potential.  LAMMPS uses a very
    old format for reading/writing EAM files.

    The potential functions are stored as they would be generated from their
    functional forms, in units eV.  However, they are written to be compliant
    with the setfl format which requires the element to be stored as V(r)*r.

    Both the electron density function and the potential energy has a 
    physical requirement that V(r_cut) = 0 and rho(r_cut) = 0.  The 
    application of a cutoff function not implemented by this class, and
    it is responsibility of the code using this class to ensure this.

    All units of measurement are in metal units.  Distance being Angstroms,
    and energy being eV.

    Args:
        None
    Attributes:
        comments(list of str):  The setfl file format has a three line comment
           section at the beginning of the file.
        filename(str):  This is the filename to either read/write too.
        symbols(list of str):  This is a list of symbols modelled by the 
           Embedded Atom Model.
        n_symbols(int): This is the number of elements for the setfl file.
        N_r(int): The number of points sampled for the distance from the nuclei
           in the electron density function, or the distance between two atoms
           in pair potentials.
        N_rho(int): The number of points where the electron density is sampled.
        d_r(float): Distance between the point where the potential functions
           and the electron density function are evaluated.
        d_rho(float): Distance between the points where the delectron density
           is sampled.
        r(numpy.ndarray):  This is a :numpy.ndarray: array of possible 
           evaluations from d_r to r_max with N_r intervals.  It does not 
           include r = 0.
        rho(numpy.ndarray): This is an :numpy.ndarray: array of possible
           evaluations from d_rho to rho_max with N_rho interals.  It does
           not include rho = 0.
        func_embedding(OrderedDict): This is a dictionary containing the 
           embedding function.  The embedding function is indexed by elements in
           the symbols attribute.  The embedding function is stored as an 
           :numpy.ndarray: of length N_r.  This attribute is initialized as None.
        func_density(OrderedDict): This is a dictionary containing the electron
           density function indexed by elements in the symbols attribute.  The
           density function is stored as an :numpy.ndarray: of length N_r.  This
           attribute is initialized as None.
        func_pair_potential(OrderedDict):  They has the format 'sym_i.sym_j' which
           are generated from unique combinations from the symbols attribute. The
           pairs are selected as s_i, s_j \in S, where i <= j.  This attribute is 
           initialized as none.
           
    Ref:
        https://sites.google.com/a/ncsu.edu/cjobrien/tutorials-and-guides/eam
        http://www.ctcms.nist.gov/potentials/
        http://lammps.sandia.gov/doc/pair_eam.html


    """

    def __init__(self):
        # some configuration attributes
        self.N_VALUES_PER_LINE_RHO = 5
        self.N_VALUES_PER_LINE_R = 5
        self.SETFL_NUM_FORMAT = "{:+24.16E}"
        self.SETFL_INT_FORMAT = "{:5d}"
        self.PAIR_KEY_FORMAT = "{}.{}"

        # these attributes should be treated as private
        self.N_headersection_lines = 5
        self.N_density_lines = None
        self.N_embedding_lines = None
        self.N_pairpotential_lines = None
        self.N_atomicsection_lines = None
        self.N_pairsection_lines = None

        # these attributes shoulde be treated as public
        self._comments = [
                "file automatically written by pypospack",
                "",
                ""]

        self.filename = None
        self._symbols = None
        self._n_symbols = None
        self.symbol_dict = None

        self.r = None
        self.rho = None
        self.func_embedding = None
        self.func_pairpotential = None
        self.func_density = None

        self._N_rho = None
        self.d_rho = None
        self.max_rho = None 
        self._N_r = None
        self.d_r = None
        self.max_r = None

        self.rcut_g = None

    @property
    def symbols(self):
        return self._symbols

    @symbols.setter
    def symbols(self,symbols):
        if isinstance(symbols,list):
            _symbols = list(symbols)
        else:
            raise ValueError("symbols must be a list of strings")

        if not all([isinstance(s,str) for s in _symbols]):
            raise ValueError("symbols must be a list of strings")

        self._symbols = _symbols
        self._n_symbols = len(self.symbols)

    @property
    def n_symbols(self):
        """
        (int): number of elements which this file describes
        """
        return self._n_symbols
    
    @property
    def N_rho(self):
        return self._N_rho

    @N_rho.setter
    def N_rho(self,nrho):
        if type(nrho) != int:
            err_msg = (
                    "N_rho must be an integer greater than zero, attempted to "
                    "set this property to {}")
            err_msg = err_msg.format(str(self.N_rho))
            raise ValueError(err_msg)
        if nrho <= 0:
            err_msg = (
                    "N_rho must be an integer greater than zero, attempted to "
                    "set this property to {}")
            err_msg = err_msg.format(str(self.N_rho))
            raise ValueError(err_msg)
        if nrho % self.N_VALUES_PER_LINE_RHO != 0:
            err_msg = (
                    "N_rho must be divisible by N_VALUES_PER_LINE_RHO, attempted "
                    "to set N_rho to {} which is not divisible by {}")
            err_msg = err_msg.format(
                    self.N_rho,
                    self.N_VALUES_PER_LINE_RHO)
            raise ValueError(err_msg)

        self._N_rho = nrho

    @property
    def N_r(self):
        return self._N_r

    @N_r.setter
    def N_r(self,nr):
        if type(nr) != int:
            err_msg = (
                    "N_r must be an integer greater than zero, attempted to "
                    "set this property to {}")
            err_msg = err_msg.format(str(nr))
            raise ValueError(err_msg)
        if nr <= 0:
            err_msg = (
                    "N_r must be an integer greater than zero, attempted to "
                    "set this property to {}")
            err_msg = err_msg.format(str(nr))
            raise ValueError(err_msg)

        if nr % self.N_VALUES_PER_LINE_R!= 0:
            err_msg = (
                    "N_r must be divisible by N_VALUES_PER_LINE_RHO, attempted "
                    "to set N_r to {} which is not divisible by {}")
            err_msg = err_msg.format(
                    nr,
                    self.N_VALUES_PER_LINE_R)
            raise ValueError(err_msg)

        self._N_r = nr
    
    @property
    def comments(self):
        return self._comments

    @comments.setter
    def comments(self,value):
        if isinstance(value,list):
            _comments = list(value)
        else:
            raise ValueError("comments must be list of length 3")
        
        if len(_comments) != 3:
            raise ValueError("comments must be a list of length 3")

        # if any comment is not a string then cast into string
        for i,c in enumerate(_comments):
            if not isinstance(c,str):
                _comments[i] = str(c)

        self._comments = list(_comments)
    
    def write(self,
            filename,
            symbols,
            r,
            rcut,
            rho,
            pair,
            embedding,
            density):
        """
        write a setfl file

        Args:
            filename(str)
            symbols(list)
            r(np.ndarray)
            rcut(numpy.ndarray)
            rho(numpy.ndarray)
            pair(collections.OrderedDict)
            embedding(collections.OrderedDict)
            density(collections.OrderedDict)
        """
        assert isinstance(filename,str)
        assert type(symbols) is list
        assert all([type(s) is str for s in symbols])
        assert isinstance(r,np.ndarray)
        assert isinstance(rcut,float)
        assert isinstance(rho,np.ndarray)
        assert type(pair) in [dict,OrderedDict]
        assert type(embedding) in [dict,OrderedDict]
        assert type(density) in [dict,OrderedDict]
        assert all([isinstance(pv,np.ndarray) for pn,pv in pair.items()])
        assert all([isinstance(ev,np.ndarray) for en,ev in embedding.items()])
        assert all([isinstance(dv,np.ndarray) for dn,ev in density.items()])
        assert all([pv.shape == r.shape for pn,pv in pair.items()])
        assert all([ev.shape == rho.shape for en,ev in embedding.items()])
        assert all([dv.shape == r.shape for dn,dv in density.items()])

        if filename is not None:
            self.filename = filename

        self.symbols = symbols
        self.r = r
        self.rho = rho
        self.pair = pair
        self.embedding = embedding
        self.density = density

        self.r_cut = rcut
        self.d_r = r[1] - r[0]
        self.N_r = r.size
        self.d_rho = rho[1] - rho[0]
        self.N_rho = rho.size

        _str_out = "".join([
            self.get_str_setfl_header_section(),
            self.get_str_setfl_atomic_section(),
            self.get_str_setfl_pairpotential_section()
            ])

        with open(self.filename,'w') as f:
            f.write(_str_out)

    def get_str_setfl_header_section(self):
        s = "".join([
            self.get_str_setfl_header_section__comments(),
            self.get_str_setfl_header_section__n_symbols_line(),
            self.get_str_setfl_header_section__nargs_line()
            ])
        return s

    def get_str_setfl_header_section__comments(self,comments=None):
        assert type(comments) in [list,type(None)]
        assert all([type(c) is str for c in comments])

        if comments is not None:
            self.comments = comments
        
        _str_out = "\n".join(self.comments)

        return str_out
    
    def get_str_setfl_header_section__n_symbols_line(self,symbols=None):
        assert type(symbols) in [list,type(None)]
        assert all([type(s) is str for s in symbols])
        
        if symbols is not None: self.symbols = symbols

        line_args = [str(self.n_symbols)] + self.symbols
        str_out = " ".join(line_args) + "\n"

        return str_out
    
    def get_str_setfl_header_section__nargs_line(self,
            N_rho=None,
            d_rho=None,
            N_r=None,
            d_r=None,
            r_cut=None):

        if N_rho is not None: self.N_rho = N_rho
        if d_rho is not None: self.d_rho = d_rho
        if N_r is not None: self.N_r = N_r
        if d_r is not None: self.d_r = d_r
        if r_cut is not None: self.r_cut = r_cut

        assert type(self.N_rho) == int
        assert type(self.d_rho) == float
        assert type(self.N_r) == int
        assert type(self.d_r) == float
        assert type(self.r_cut) == float

        #"{0:5d}{1:24.16e}{2:5d}{3:24.16e}{4:24.16f}"
        str_out = "{}{}{}{}{}\n".format(
                self.SETFL_INT_FORMAT,self.SETFL_NUM_FORMAT,
                self.SETFL_INT_FORMAT,self.SETFL_NUM_FORMAT,
                self.SETFL_NUM_FORMAT).format(
                    self.N_rho,self.d_rho,
                    self.N_r,self.d_r,
                    self.r_cut)

        return str_out

    def get_str_setfl_atomic_section(self):
        _list_str = []
        for s in self.symbols:
            _list_str.append(self.get_str_setfl__atomic_description(s))
            _list_str.append(self.get_str_setfl__embeddding_function(s))
            _list_str.append(self.get_str_setfl__density_function(s))
        _str_out = "\n".join(_list_str)
        return _str_out

    def get_str_setfl__atomic_descriptions(self,symbol):
        _atomic_information = OrderedDict()        
        if symbol == 'Ni':
            _atomic_information['an'] = 28   
            _atomic_information['amu'] = 5.871
            _atomic_information['a0'] =  3.518121
            _atomic_information['latt_type'] = 'fcc'
        else:
            raise ValueError("symbol not in database")
        
        str_out = "".format([
                self.SETFL_INT_FORMAT.format(_atomic_information['an']),
                self.SETFL_NUM_FORMAT.format(_atomic_information['amu']),
                self.SETFL_NUM_FORMAT.format(_atomic_information['a0']),
                " " + _atomic_information['fcc'],"\n"
                ])
    
    def get_str_setfl__embedding_function(self,symbol):

        _embed = self.embedding[s]
        _sz = _embed.size
        _n = self.N_VALUES_PER_LINE_RHO
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_embed[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)

        return _str_out

    def get_str_setfl__density_function(self,symbol):
        
        _dens = self.density[s]
        _sz = _dens.size
        _n = self.N_VALUES_PER_LINE_R
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_dens[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)

        return _str_out

    def get_str_setfl_pairpotential_section(self):
        _str_out = ""
        for i1,s1 in self.symbols:
            for i2,s2 in self.symbols:
                if i1 <= i2:
                    _str_out += get_str_setfl__pair_function(s1,s2)
        return _str_out

    def get_str_setfl__pair_function(self,symbol1,symbol2):
        pn = "{}{}".format(symbol1,symbol2)
        _pair = self.pair[pn] * self.r
        _sz = _dens.size
        _n = self.N_VALUES_PER_LINE_R
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_pair[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)
    
    def read(self,
            filename=None,
            autoprocess=True):
        if filename is not None:
            assert type(filename) is str
            self.filename = filename

        self.lines = filesystem.read_file_as_lines(self.filename)

        assert isinstance(self.lines,list)
        if autoprocess is True:
            self.process_setfl_header_section()
            self.process_setfl_atomic_section()
            self.process_setfl_pairpotential_section()

    def process_setfl_header_section(self):
        """
        Process the header section of the setfl file.

        The header section consists of two sections
        """
        self.process_setfl_comments()
        self.process_setfl_n_symbols()
        self.process_setfl_rho_r_args()

        # determine the rho and r arrays
        self.max_rho = self.N_rho * self.d_rho
        self.max_r = self.N_r * self.d_r
        if self.max_r <= self.r_cut:
            self.r_cut = self.max_r

        self.r = self.max_r*np.linspace(1,self.N_r,self.N_r)/self.N_r
        self.rho = self.max_rho*np.linspace(1,self.N_rho,self.N_rho)/self.N_rho
   
        if self.N_rho % self.N_VALUES_PER_LINE_RHO != 0:
            raise ValueError("N_rho not divisible by {}".format(
                self.N_VALUES_PER_LINE_RHO))
        if self.N_r % self.N_VALUES_PER_LINE_R != 0:
            raise ValueError("N_r not divisible by {}".format(
                self.N_VALUES_PER_LINE_R))

        N_density_lines = self.N_rho / self.N_VALUES_PER_LINE_RHO 
        N_embedding_lines = self.N_r / self.N_VALUES_PER_LINE_R
        N_pairpotential_lines = self.N_r / self.N_VALUES_PER_LINE_R

        # needs to be cast as a integer
        N_density_lines = int(N_density_lines)
        N_embedding_lines = int(N_embedding_lines)
        N_pairpotential_lines = int(N_pairpotential_lines)

        # set these values as attributes
        self.N_density_lines = N_density_lines
        self.N_embedding_lines = N_embedding_lines
        self.N_pairpotential_lines = N_embedding_lines

        self.N_atomicsection_lines = self.n_symbols*(
                1+self.N_density_lines+self.N_embedding_lines)

        pairs = []
        for i1,s1 in enumerate(self.symbols):
            for i2,s2 in enumerate(self.symbols):
                if i1 >= i2:
                    pairs.append([s1,s2])

        self.N_pairs = len(pairs)
        self.N_pairsection_lines = self.N_pairs * N_pairpotential_lines 

    def process_setfl_comments(self):
        """

        The header section of the setfl files containes three lines which
        can have arbitrary text in it.

        Args:
            None
        Returns:
            (list of str)
        """
        comment_lines = [0,1,2]
        comments = []
        for n in comment_lines:
            comments.append(self.lines[n])
        self.comments = list(comments)

        return comments

    def process_setfl_n_symbols(self):
        line = self.lines[3]
        args = line.split()
        
        n_symbols = int(args[0])
        symbols = []
        for i in range(1,len(args)):
            symbols.append(args[i])

        if n_symbols != len(symbols):
            err_msg = (
                    "The setfl file indicates a different number of atoms than "
                    "number of species provided")
            # -- debugging information
            print("line:\n",line)
            print("symbols:",symbols)
            print("len(symbols):",len(symbols))
            print("n_symbols:",n_symbols)
            raise ValueError(err_msg)

        try:
            self.symbols = symbols
        except ValueError as e:
            # -- debugging information
            print("line:\n",line)
            print("symbols:",symbols)
            print("len(symbols):",len(symbols))
            print("n_symbols:",n_symbols)
            raise

        return list(self.symbols), self.n_symbols

    def process_setfl_rho_r_args(self,string=None):
        # initialize private variables
        line = args = None

        # process method arguments
        if type(string) is str:
            line = string
        else:
            line = self.lines[4]

        # split the line
        args = line.split()

        self.N_rho = int(args[0])
        self.d_rho = float(args[1])
        self.N_r = int(args[2])
        self.d_r = float(args[3])
        self.r_cut = float(args[4])

        return_dict = {
                'N_rho':self.N_rho,
                'd_rho':self.d_rho,
                'N_r':self.N_r,
                'd_r':self.d_r,
                'r_cut':self.r_cut}

        return return_dict

    def process_setfl_atomic_section(self):
        # initialize
        self.symbol_dict = OrderedDict()
        self.func_embedding = OrderedDict()
        self.func_density = OrderedDict()

        # Local copy of attributes
        start_line = self.N_headersection_lines + 1
        N_embedding_lines = self.N_embedding_lines
        N_density_lines = self.N_density_lines
        N_atomicsection_lines = self.N_atomicsection_lines
        N_atomicdefinition_lines = 1
        for i_sym, sym in enumerate(self.symbols):
            # NOTICE: The use of the notion self.lines[i-1] 
            # I decided to calculate the line we are on using an index of one,
            # then using (i-1) to convert to the python convertion of an
            # index starting a 0

            # read atomic_definition
            i = start_line + i_sym*(1+N_density_lines+N_embedding_lines)
            args = self.lines[i-1].split()
            self.symbol_dict[sym] = OrderedDict({
                    'atomic_number':int(args[0]),
                    'atomic_mass':float(args[1]),
                    'lattice_constant':float(args[2]),
                    'lattice_type':args[3]})
            
            # read embedding function
            f_embed = []
            for i_line in range(N_density_lines):
                i = start_line \
                        + i_sym*N_atomicsection_lines\
                        + N_atomicdefinition_lines \
                        + i_line
                args = self.lines[i-1].strip().split()
                f_embed += [float(arg) for arg in args]
                self.func_embedding[sym] = np.array(f_embed)
            
            # read density function
            f_dens  = []
            for i_line in range(N_embedding_lines):
                i = start_line \
                        + i_sym*N_atomicsection_lines\
                        + N_atomicdefinition_lines\
                        + N_embedding_lines\
                        + i_line
                args = self.lines[i-1].strip().split()
                f_dens += [float(arg) for arg in args]
            self.func_density[sym] = np.array(f_dens)

    def process_setfl_pairpotential_section(self):
        """
        
        The SETFL convention for storing V(r) is as r*V(r).  This was probably 
        important during the time when computational cycles were expensive.  
        However, we divide this storage value by r to get the potential in 
        eV.
        """
        start_line = self.N_headersection_lines + self.N_atomicsection_lines + 1
        N_pair_lines = self.N_pairpotential_lines

        pairs = []
        for idx1,sym1 in enumerate(self.symbols):
            for idx2,sym2 in enumerate(self.symbols):
                if idx1 <= idx2:
                    pairs.append([sym1,sym2])

        self.func_pairpotential = OrderedDict()
        for i_pair,pair in enumerate(pairs):

            f_pair = []
            for i_line in range(N_pair_lines):
                i = start_line + i_pair * N_pair_lines + i_line
                args = self.lines[i-1].strip().split()
                f_pair += [float(arg) for arg in args]
            k = self.PAIR_KEY_FORMAT.format(pair[0],pair[1])
            self.func_pairpotential[k] = np.array(f_pair)/self.r


    def calc_N_lines_atomic_section(self,
            N_symbols=None,
            N_rho=None,
            N_r=None):

        # declare local variables
        _N_symbols = N_symbols
        _N_rho = N_rho
        _N_r = N_r

        # process argument N_symbols
        if _N_symbols is None: _N_symbols = self.n_symbols
        if _N_rho is None: _N_rho = self.N_rho
        if _N_r is None: _N_r = self.N_r

        # check these assertions
        # TODO: recode these to raise exceptions if errors pop up here
        assert type(_N_rho) is int
        assert _N_rho > 0
        assert _N_r % self.N_VALUES_PER_LINE_R == 0
        assert type(_N_r) is int
        assert _N_rho > 0
        assert _N_rho % self.N_VALUES_PER_LINE_RHO == 0

        N_lines_atomic_section = \
                N_symbols * (
                        self.N_rho / self.N_LINES_PER_LINE_RHO
                        + self.N_r / self.N_LINES_PER_LINE_R)

        return N_lines_atomic_section

    def calc_N_lines_potential_section(self,
            N_symbols=None,
            N_r=None):

        # declare some local variables
        _N_symbols=N_symbols
        _N_r=N_r

        # process arguments
        if _N_symbols is None: _N_symbols = self.n_symbols
        if _N_r is None: _N_r = self.N_r

        assert type(_N_symbols) is int
        assert _N_symbols > 0
        assert type(_N_r) is int
        assert _N_r > 0
        assert _N_r % self.N_LINES_PER_LINE_R == 0
   
        N_pairs = 0
        for i1 in range(self.n_symbols):
            for j1 in range(self.n_symbols):
                N_pairs += 1

        N_lines_potential_section = \
                N_pairs * (self.N_r /self.N_LINES_PER_LINE_R)

        return N_lines_potential_section

    def get_setfl_head(self):
        s  = self.get_comments_as_string()
        s += self.get_n_symbols_line_as_string()
        s += self.get_setfl_nargs_as_string()
        return s


    def get_n_symbols_line_as_string(self,symbols=None):
        if symbols is not None:
            self.symbols = symbols

        if self._symbols is None:
            msg_out = "the symbols attribute has not been initialized"
            raise ValueError(msg_out)

        line_args = [str(self.n_symbols)] + self._symbols
        str_out = " ".join(line_args) + "\n"

        return str_out

    def get_setfl_nargs_line_as_string(self,
            N_rho=None,
            d_rho=None,
            N_r=None,
            d_r=None,
            r_cut=None):

        if N_rho is not None: self.N_rho = N_rho
        if d_rho is not None: self.d_rho = d_rho
        if N_r is not None: self.N_r = N_r
        if d_r is not None: self.d_r = d_r
        if r_cut is not None: self.r_cut = r_cut

        assert type(self.N_rho) == int
        assert type(self.d_rho) == float
        assert type(self.N_r) == int
        assert type(self.d_r) == float
        assert type(self.r_cut) == float

        #"{0:5d}{1:24.16e}{2:5d}{3:24.16e}{4:24.16f}"
        str_out = "{}{}{}{}{}\n".format(
                self.SETFL_INT_FORMAT,self.SETFL_NUM_FORMAT,
                self.SETFL_INT_FORMAT,self.SETFL_NUM_FORMAT,
                self.SETFL_NUM_FORMAT).format(
                    self.N_rho,self.d_rho,
                    self.N_r,self.d_r,
                    self.r_cut)

        return str_out
        
    def get_setfl_func_density_to_string(self,symbol):
        N = self.N_VALUES_PER_LINE_RHO
        dens = self.func_density[symbol].tolist()
        
        str_out = ""
        for i in range(len(dens)/N + 1):
            N_start = i*N
            N_stop = (i+1)*N
            str_out += self.SETFL_NUM_FORMAT.join(
                    dens[N_start:N_stop]) + "\n"
        return str_out

    def get_setfl_func_embedding_to_string(self,symbol):
        N = self.N_VALUES_PER_LINE_R
        embed = self.func_density[symbol].tolist()

        str_out = ""
        for i in range(len(embed)/N + 1):
            N_start = i*N
            N_stop = (i+1)*N
            str_out += self.SETFL_NUM_FORMAT.join(
                    embed[N_start:N_stop]) + "\n"
        return str_out

    def get_setfl_density_function_to_string(self,symbol1,symbol2):
        N = self.N_VALUES_PER_LINE_R
        key = self.PAIR_KEY_FORMAT.format(symbol1,symbol2)
        pair = self.func_pair_potential[k].tolist()

        str_out = ""
        for i in range(len(pair)/N + 1):
            N_start = i*N
            N_stop = (i+1)*N
            str_out += self.SETFL_NUM_FORMAT.join(
                    pair[N_start:N_stop]) + "\n"
        return str_out

class EamFile2(object):
    """

    Ref:
        https://sites.google.com/a/ncsu.edu/cjobrien/tutorials-and-guides/eam
        http://www.ctcms.nist.gov/potentials/
        http://lammps.sandia.gov/doc/pair_eam.html
    """
    def __init__(self):
        self._comments = []
        self._comments.append("file automatically written by pyposmat")
        self._comments.append("")
        self._comments.append("")

        self._title = None
        self._elist = None
        self._elements = None
        self._fname = None

        self._r = None
        self._embed = {}
        self._pair = {}
        self._dens = {}

        self._Nrho = None
        self._drho = None
        self._Nr = 500
        self._dr = 500

        self._rcut_g = None

    @property
    def comments(self):
        return self._comments

    @comments.setter
    def comments(self, comments):
        if len(comments) <= 3:
            for i,v in enumerate(comments):
                self._comments[i] = v.strip()
        else:
            raise Exception("comments must have less than three items")

    @property
    def title(self):
        return self._title

    @property
    def element_list(self):
        return self._elist

    @property
    def filename(self):
        return self._fname

    @property
    def embed(self):
        return self._embed

    @property
    def pair(self):
        return self._pair

    @property
    def dens(self):
        return self._dens

    @property
    def Nrho(self):
        return self._Nrho

    @Nrho.setter
    def Nrho(self, Nrho):
        self._Nrho = Nrho

    @property
    def drho(self):
        return self._drho

    @drho.setter
    def drho(self, drho):
        self._drho = drho

    @property
    def Nr(self):
        return self._Nr

    @Nr.setter
    def Nr(self, Nr):
        self._Nr = Nr

    @property
    def rcut_global(self):
        return self._rcut_g

    @rcut_global.setter
    def rcut_global(self, rcut_g):
        self._rcut_g = rcut_g

    def write(self,filename,filetype = 'setfl'):
        if filetype == 'setfl':
            self._write_setfl()

    def _write_setfl(self):
        self.file = open(self._fname, 'w')
        self._write_eam_head()
        self._write_singleeam()
        self._write_paiream()
        self.file.close()

    def _write_eam_head(self):
        for v in self._comments:
            self.file.write(v+"\n")

        # write number of atoms
        n_elements = len(self.elist)
        line = '    {}'.format(n_elements)
        for s in self.elist:
            line += (' {0:s}'.format(s))
        line += '\n'
        self.file.write(line)
        
        # 
        self.file.write('{0:5d}{1:24.16e}{2:5d}{3:24.16e}{4:24.16f}\n'.format(\
                             self.Nrho, self.drho,
                             self.Nr,   self.dr,
                             self.globalcutoff))

    def _write_atomic_section(self):
        for el in self.elist:
            an = toe.get_atomic_number(el)
            am = toe.get_atomic_mass(el)
            lc = toe.get_lattice_constant(el)
            lt = toe.get_lattice_type(el)
            line = ("{}{}{}{}\n".format(an,am,lc,lt))
            self._write_embed_func('el')
            self._write_dens_func('el')
    
    def _write_potential_section(self):
        N_el = len(self.elist)
        for i_el in range(N_el):
            for j_el in range(N_el):
                if i_el <= j_el:
                    self._write_potential_func('{}{}'.format(self._elist[i_el],
                                                             self._elist[j_el]))

class EamFileWriter(object):
    # Uses functions developed by O'Brien and Foiles
    def __init__(self, title, elements, filename, Nrho,drho,Nr,dr):
        self.comment1 = "file automatically written by pyposmat\n"
        self.comment2 = "\n"
        self.comment3 = "\n"
        self.title = title
        self.elist = elements
        self.elements = Element(elements[0])
        self.filename     = filename
        
        # EAM functions
        self.fembed = []
        self.frho   = []
        self.fpairr = []

        self.Nrho = Nrho
        self.drho = drho
        self.Nr   = Nr
        self.dr   = dr

        self.globalcutoff = self.Nr * self.dr
        
    def write(self, file_type = 'setfl'):
        if file_type == 'setfl':
            self._write_setfl()
            
    def compute_embedding(self, func_embedding, params):
        self.rho = np.linspace(start = 0, 
                          stop = self.Nrho * self.drho, 
                          num = self.Nrho)
        if func_embedding.__name__ == 'func_embedding_universal':
            self.elements.embed =  func_embedding(self.rho, *params)
        else:
            raise Exception("{} is not an acceptable embedding function".format(func_embedding.__name__))

    def compute_density(self, func_density, params):
        self.r   = np.linspace(start = 0, 
                          stop = self.Nr * self.dr, 
                          num = self.Nr)
        if func_density.__name__ == "func_edensity_3s_metal":
            self.elements.rho =  func_density(self.r, *params)
        else:
            raise Exception("{} is not an acceptable density function".format(func_density.__name__))
            
    def compute_pairr(self, func_pairr, params):
        self.r   = np.linspace(start = 0, 
                          stop = self.Nr * self.dr, 
                          num = self.Nr)
        if func_pairr.__name__ == "func_pairr_effective_charge":
            self.pairr = func_pairr(self.r, *params)
        else:
            raise Exception("{} is not an acceptable pair function".format(func_pairr.__name__))
            
    def _write_setfl(self):
        self.file = open(self.filename, 'w')
        self.write_eamhead()
        self.write_singleeam()
        self.write_paiream()
        self.file.close()

    def write_eamhead(self):
        ''' write header for eam file'''
        self.file.write(self.comment1)
        self.file.write(self.comment2)
        self.file.write(self.comment3)
        #write num atoms
        if len(self.elist)==1:
            self.file.write('    1 {0:s}\n'.format(self.elist[0]))
        elif len(self.elist)==2:
            self.file.write('    2 {0:s} {1:s}\n'.format(self.elist[0],self.elist[1]))
        else:
            pass
            #raise(StandardError, "Only two elements can be used at this time!!!")

        self.file.write('{0:5d}{1:24.16e} {2:5d}{3:24.16e}{4:24.16f}\n'.format(\
                             self.Nrho, self.drho,
                             self.Nr,   self.dr,
                             self.globalcutoff))

    def write_singleeam(self):
        ''' Write Single atom section of  EAM/ALLOY style (funcfl) Potential for LAMMPS'''

        self.file.write('{0:5d}{1:15.5g}{2:15.5g} {3:s}\n'.format(self.elements.atnum,
                                                            self.elements.mass,
                                                            self.elements.alat,
                                                            self.elements.lattype))
        #write embedding fxn for Nrho values, with 5 values on each row
        n_val_per_row = 5
        i = 0
        while i < self.Nrho:                      #<---- track rho[i] per TOTAL
            j = 0
            str_out = ""
            while j < n_val_per_row:              #<---- track rho[j] per ROW
                str_out += "{:24.16e}".format(self.elements.embed[i])
                i += 1
                j += 1
            str_out += "\n"
            self.file.write(str_out)
            
        #write density for Nr values, with 5 values on each row
        n_val_per_row = 5
        i = 0
        while i < self.Nr:
            j = 0
            str_out = ""
            while j < n_val_per_row:
                str_out += "{:24.16e}".format(self.elements.rho[i])
                i += 1
                j += 1
            str_out += "\n"
            self.file.write(str_out)
                

    def write_paiream(self):
        ''' Write pair section of EAM/ALLOY style (funcfl) Potential for LAMMPS'''
        n_val_per_row = 5
        i = 0
        while i < self.Nr:
            j = 0
            str_out = ""
            while j < n_val_per_row:
                str_out += "{:24.16e}".format(self.pairr[i])
                i += 1
                j += 1
            str_out += "\n"
            self.file.write(str_out)

class EamFileReader:
    def __init__(self, filename):
      self.filename = filename

    def read(self):
        if self.filename.endswith('eam'):
            self._read_funcfl_file()
        else: 
            self._read_setfl_file()

    def _read_funcfl_file(self):
        rphi = np.sqrt(27.2*0.529)
        eV2J = spc.value('electron volt')

        print('Reading EAM potential in single-element (funcfl) format')
        f = open(self.filename, 'r')

        line = f.readline()  # comment line
        self.atomic_number = int(f.readline().split()[0])

        line = f.readline().split()
        self.nrho    =   int(line[0]) 
        self.drho    = float(line[1])
        self.nr      =   int(line[2])
        self.dr      = float(line[3])
        self.cutoff  = float(line[4])

        self.fembed  = np.zeros((self.nrho,2))
        self.effchg  = np.zeros((self.nr,2))
        self.fdens   = np.zeros((self.nr,2))

        # read embedding functional
        line    = f.readline().split()
        n_col   = len(line)
        n_lines = int(self.nrho/self.ncol)-1

        i_rows = 0
        while i_rows < self.nrho:
            for val in line:
                self.fembed[i_rows,0] = i_rows*self.drho
                self.fembed[i_rows,1] = float(val)
                i_rows += 1
                if i_rows >= self.nrho: 
                    break
            line=f.readline().split()
        
        # read charge function
        i_rows = 0
        while i_rows < self.nr:
            for val in line:
                self.effchg[i_rows,0]=i_rows*self.dr
                self.effchg[i_rows,1]=float(val)*rphi/(1.e-30+self.effchg[self.i_rows,0])
                i_rows += 1
                if i_rows >= self.nr: 
                    break
            line=f.readline().split()

        # read density function
        i_rows = 0
        while i_rows < self.nr:
            for val in line:
                self.fdens[i_rows,0]=i_rows*self.dr
                self.fdens[i_rows,1]=float(val)
                i_rows += 1
                if i_rows >= self.nr:
                    break
            line=f.readline().split()

    def _read_setfl_file(self):
        print('Reading EAM potential in EAM/ALLOY miltielement format (setfl)')
        f=open(self.filename,'r')
        self.comment_1 = f.readline()
        self.comment_2 = f.readline()
        self.comment_3 = f.readline()
        
        line = f.readline().split() #line 4
        self.n_elements  = int(line[0])
        self.element_list = line[1:]

        line = f.readline().split() #line 5
        self.Nrho =     int(line[0]) 
        self.drho =   float(line[1])
        self.Nr =       int(line[2])
        self.dr =     float(line[3])
        self.cutoff = float(line[4])

        self.fembed = {}
        self.fdens  = {}
        self.effpot = {}
        
        # initialize embedding potentials
        for element_i in self.element_list:
            self.fembed[element_i] = np.zeros((self.Nrho, 2))
            
        # initialize electron densities
        for element_i in self.element_list:
            self.fdens[element_i]  = np.zeros((self.Nr,   2))

        # initialize pair potentials
        for element_i in self.element_list:
            for element_j in self.element_list:
                pot_name = "{}{}".format(element_j,element_i)
                if not(pot_name in self.effpot):
                    pot_name = "{}{}".format(element_i,element_j)
                    self.effpot[pot_name] = np.zeros((self.Nr,2))
                
        #read embedding function and density function for each element
        for element in self.element_list:

            # read comment line
            line = f.readline().split()
            atomic_number    =   int(line[0]) #line6
            atomic_mass      = float(line[1])
            lattice_constant = float(line[2])
            lattice_type     =       line[3]
            
            #read embedding functional
            #print(line)
            #ncol    = len(line)
            #n_lines = int(self.nrho/ncol)-1

            print("Reading embedding energy of element: {}".format(element))                
            i_val = 0 #tracks number of values of rho
            line    = f.readline().split()
            while i_val < self.Nrho:
                for val in line:
                    self.fembed[element][i_val,0] = i_val*self.drho
                    self.fembed[element][i_val,1] = float(val)
                    i_val += 1
                line=f.readline().split()

            print("Reading density function of element: {}".format(element))                        
            #read density function
            i_val = 0 #track number of rows of density fxn
            while i_val < self.Nr:
                for val in line:
                    self.fdens[element][i_val,0]=i_val*self.dr
                    self.fdens[element][i_val,1]=float(val)
                    i_val+=1
                if i_val==self.Nr:
                    continue
                else:
                    line=f.readline().split()
    
        #read pair potential functions
        for element_i in self.element_list:
            for element_j in self.element_list:
                pot_name = "{}{}".format(element_i,element_j)
                if pot_name in self.effpot:
                    print("Reading pair potential: {}".format(pot_name))
                    i_val = 0 
                    while i_val < self.Nr:
                        for val in line:
                            self.effpot[pot_name][i_val,0] = float(i_val)*self.dr
                            self.effpot[pot_name][i_val,1] = float(val)
                            i_val += 1
                        line=f.readline().split()
                        self.effpot[pot_name][0,1:]=0.

class Element:

  def __init__(self,symbol):
    self.setPropertiesBySymbol(symbol)

  def setProperties(self,an,n,sym,m,lt,a,co,embedvals,rhovals,potvals):
    self.atnum = an
    self.name  = n
    self.symbol = sym
    self.mass = m
    self.lattype = lt
    self.alat = a
    self.cutoff = co
    self.embed = embedvals
    self.rho = rhovals

  def setPropertiesBySymbol(self,symbol):
    properties = ''
    if symbol == "Ni":
      properties = (28,"Nickel","Ni",58.6934,'FCC',3.52,0,[],[],[])
    self.setProperties(*properties)
