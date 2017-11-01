import numpy as np
import scipy.constants as spc
# import pyflamestk.potentials
#import pyflamestk.tableofelements as toe
class EamCurveFitter:
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

class EamSetflFile(object):
    """

    Writes an Setfl file
    Ref:
        https://sites.google.com/a/ncsu.edu/cjobrien/tutorials-and-guides/eam
        http://www.ctcms.nist.gov/potentials/
        http://lammps.sandia.gov/doc/pair_eam.html
    """

    def __init__(self):
        self.comments = [
                "file automatically written by pypospack",
                "",
                ""]

        self.title = None
        self.elist is None
        self.elements is None

        self.r = None
        self.embed = None
        self.pair = None
        self.dens = None

        self.N_rho = None
        self.d_rho = None
        self.Nr = 500
        self.dr = 500

        self.rcut_g = None
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

    def write(self,filetype = 'setfl'):
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
