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
                "--- this line empty ---",
                "--- this line empty ---"]

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
        assert all([isinstance(dv,np.ndarray) for dn,dv in density.items()])
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

        _str_out = "\n".join([
            self.get_str_setfl_header_section().strip(),
            self.get_str_setfl_atomic_section().strip(),
            self.get_str_setfl_pairpotential_section().strip()
            ])

        with open(self.filename,'w') as f:
            f.write(_str_out)

    def get_str_setfl_header_section(self):
        _str_out = "\n".join([
            self.get_str_setfl_header_section__comments(),
            self.get_str_setfl_header_section__n_symbols_line(),
            self.get_str_setfl_header_section__nargs_line()
            ])
        return _str_out

    def get_str_setfl_header_section__comments(self,comments=None):
        assert type(comments) in [list,type(None)]

        if comments is not None:
            self.comments = comments
        assert all([type(c) is str for c in self.comments])
        _str_out = "\n".join(self.comments)

        return _str_out
    
    def get_str_setfl_header_section__n_symbols_line(self,symbols=None):
        assert type(symbols) in [list,type(None)]
        
        if symbols is not None: self.symbols = symbols
        assert all([type(s) is str for s in self.symbols])

        line_args = [str(self.n_symbols)] + self.symbols
        _str_out = " ".join(line_args)

        return _str_out
    
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
        assert isinstance(self.d_rho,float)
        assert type(self.N_r) == int
        assert isinstance(self.d_r,float)
        assert isinstance(self.r_cut, float)

        #"{0:5d}{1:24.16e}{2:5d}{3:24.16e}{4:24.16f}"
        _str_out = " ".join([
                self.SETFL_INT_FORMAT.format(self.N_rho),
                self.SETFL_NUM_FORMAT.format(self.d_rho),
                self.SETFL_INT_FORMAT.format(self.N_r),
                self.SETFL_NUM_FORMAT.format(self.d_r),
                self.SETFL_NUM_FORMAT.format(self.r_cut)])

        return _str_out

    def get_str_setfl_atomic_section(self):
        _list_str = []
        for s in self.symbols:
            _list_str.append(self.get_str_setfl__atomic_description(s))
            _list_str.append(self.get_str_setfl__embedding_function(s)[:-1])
            _list_str.append(self.get_str_setfl__density_function(s)[:-1])
        _str_out = "\n".join(_list_str)
        return _str_out

    def get_str_setfl__atomic_description(self,symbol):
        _atomic_information = OrderedDict()        
        if symbol == 'Ni':
            _atomic_information['an'] = 28   
            _atomic_information['amu'] = 5.871
            _atomic_information['a0'] =  3.518121
            _atomic_information['latt_type'] = 'fcc'
        elif symbol == 'Al':
            _atomic_information['an'] = 13
            _atomic_information['amu'] = 2.698
            _atomic_information['a0'] = 4.046
            _atomic_information['latt_type'] = 'fcc'
        else:
            raise ValueError("symbol not in database")
        
        _str_out = "".join([
                self.SETFL_INT_FORMAT.format(_atomic_information['an']),
                self.SETFL_NUM_FORMAT.format(_atomic_information['amu']),
                self.SETFL_NUM_FORMAT.format(_atomic_information['a0']),
                " {:s}".format(_atomic_information['latt_type'])
                ])
        return _str_out

    def get_str_setfl__embedding_function(self,symbol):

        _embed = self.embedding[symbol]
        _sz = _embed.size
        _n = self.N_VALUES_PER_LINE_RHO
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_embed[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)

        return _str_out

    def get_str_setfl__density_function(self,symbol):
        
        _dens = self.density[symbol]
        _sz = _dens.size
        _n = self.N_VALUES_PER_LINE_R
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_dens[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)

        return _str_out

    def get_str_setfl_pairpotential_section(self):
        _str_out = ""
        for i1,s1 in enumerate(self.symbols):
            for i2,s2 in enumerate(self.symbols):
                if i1 <= i2:
                    _str_out += self.get_str_setfl__pair_function(s1,s2)
        return _str_out

    def get_str_setfl__pair_function(self,symbol1,symbol2):
        pn = "{}{}".format(symbol1,symbol2)
        _pair = self.pair[pn] * self.r
        _sz = _pair.size
        _n = self.N_VALUES_PER_LINE_R
        
        _format = _n*self.SETFL_NUM_FORMAT + "\n"
        _lines = [_format.format(*_pair[i:i+5]) for i in range(0,_sz,_n)]
        _str_out = "".join(_lines)
   
        return _str_out

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

