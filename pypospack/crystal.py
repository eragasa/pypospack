# -*- coding: utf-8 -*-
"""Crystallographic classes and methods for pypospack.

This module represents simulation cells are crystallographic representions and 
classes and methods to support this represention.

Attributes:
    iso_chemical_symbols(list): list of chemical symbols
    atom_info(dict): dictionary of chemical symbols, and a dictionary of their values
"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import copy, subprocess, yaml, os
import numpy as np
import numpy.linalg as linalg
import ase.atoms

iso_chem_symbols = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
        'Si','P','S','Cl','Ar','K', 'Ca','Sc','Ti','V', 'Cr','Mn','Fe','Co',
        'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr','Nb',
        'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs',
        'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
        'Yb','Lu','Hf','Ta','W', 'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi',
        'Po','At','Rn','Fr','Ra','Ac','Th','Pa','U', 'Np','Pu','Am','Cm','Bk',
        'Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg',
        'Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo']

atom_info = {}
atom_info['Ar'] = {'atomic_mass':39.948}
atom_info['Mg'] = {'atomic_mass':24.305}
atom_info['O'] = {'atomic_mass':15.9994}
atom_info['Si'] = {'atomic_mass':28.0855}
atom_info['Ni'] = {'atomic_mass':58.6934}
atom_info['Al'] = {'atomic_mass':26.981539}

def get_amu(symbol):
    """ get atomic mass unit

    this method gets the the atomic mass unit for a given chemical symbol.

    Args:
        symbol(str): the ISO chemical symbol of an atom

    Returns:
        float: the mass of the atom, in atomic mass units
    """

    return atom_info[symbol]['atomic_mass']

def cartesian2direct(x,H):
    """ transforms cartesian coordinates to direct coordinates

    Args:
        x (numpy.ndarray): the cartesian coordinates
        H (numpy.ndarray): the H matrix of lattice vectors

    Return:
        (numpy.ndarray): the direct coordinate vector
    """
    if isinstance(x,list):
        x = np.array(x)
        x.shape = (3,1)
    elif isinstance(x,np.ndarray):
        x = np.copy(x)
        x.shape = (3,1)

    u = np.dot(linalg.inv(H.T),x)
    u.shape = (3,)
    return u

class StructureDatabase(object):
    """ structure database 

    Attributes:
        filename(str): file to read/write the yaml file
        directory(str): the directory of the structure database
        structures(dict): key is structure name, value is another dict with 
            key/value pairs for 'filename' and 'filetype'
    """

    def __init__(self):
        self.filename = 'pypospack.structure.yaml'
        self.directory = None
        self.structures = {}

    def add_structure(self,name,filename,filetype):
        self.structures[name] = {'filename':filename,'filetype':filetype}

    def contains(self,structure):
        """ check to see that the structure in the structure database

        Args:
            structure(str):

        Returns:
            bool: True if structure is in structure database. False if the 
                structure is not in the structure database
        """

        return structure in self.structures.keys()

    def read(self,fname = None):
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed 
                then use the filename attribute.  If the filename is set, then 
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            self.structure_db = yaml.load(open(self.filename))
        except:
            raise

        self.directory = self.structure_db['directory']
        self.structures = copy.deepcopy(self.structure_db['structures'])

    def write(self,fname = None):
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.structure_db = {}
        self.structure_db['directory'] = self.directory
        self.structure_db['structures'] = {}
        self.structure_db['structures'] = copy.deepcopy(self.structures)

        # dump to as yaml file
        with open(fname,'w') as f:
            yaml.dump(self.structure_db, f, default_flow_style=False)

    def check(self):
        """sanity checks for the fitting database
        
        This method checks the fitting database for the following errors: (1)
        the structure database exists, (2) files for the structure database
        exist

        Returns:
            str: returns a string if there is a problem
        Raises:
            ValueError: if there is a problem with the configuration
        """

        src_dir = self.directory
       
        # check to see if the source directory exists
        if not os.path.exists(src_dir):
            err_msg = "cannot find simulation directory\n"
            err_msg += "\tcurrent_working_directory:{}\n".format(os.getcwd())
            err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
            return err_msg
        
        # check to see if the source directory is a directory
        if not os.path.isdir(src_dir):
            err_msg = "path exists, is not a directory\n"
            err_msg += "\tcurrent_working_directory:{}".format(os.getcwd())
            err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
            return err_msg

        # check to see if files exist in the source directory
        files_exist = True
        msg = "structure files are missing:\n"
        for name, v in self.structures.items():
            filename = os.path.join(src_dir,v['filename'])
            if not os.path.isfile(filename):
                files_exist = False
                msg += "\t{}:{}\n".format(name,filename)

        if not files_exist:
            return msg
        else:
            return True

    def get_structure_dict(self,name):
        """

        Args:
            name(str): name of the structure
        """

        structure_db_dir = self.directory
        structure_filename = self.structures[name]['filename']
        structure_dict = {}
        structure_dict['name'] = name
        structure_dict['filename'] = os.path.join(
                structure_db_dir,
                structure_filename)

        return copy.deepcopy(structure_dict)

class Atom(object):
    """description of an atom

    This position is a data structure which contains information about an
    individual atom
    
    Args:
        symbol (str): the standard ISO symbol for an element
        position (list of float): the position of the atom the units
           are dependent upon use.
        magmom (float): the magnetic moment of the atom.
    
    Attributes:
        symbol (str): the standard ISO symbol for an element
        position (numpy.ndarray): the position of the atom, usually in direct
            coordinates
        magentic_moment (float): the magnetic moment of the atom
    """
    def __init__(self, symbol, position, magmom = 0):

        self.symbol = symbol
        if isinstance(position,list):
            self.position = np.array(position)
        elif isinstance(position,np.ndarray):
            self.position = position.copy()
        else:
            raise TypeError('position must either be a list of numeric values or a numpy array')
        self.magnetic_moment = magmom
       

class SimulationCell(object):
    """A structural representation of a material system
  
    A structural system consists of a vector of lattice vector which forms the 
    boundaries of the simulation cell, and an atomic basis of atoms defined in
    the direct coordinates of the lattice vectors.

    Args:
        obj (optional): if this argument is set then the this constructor acts 
            as a copy constructor.  Will takse :obj:`ase.atoms.Atoms` and

    Attributes:
        comment (str): a descriptive description of the crystal structure
        atomic_basis (:obj:`list` of :obj:`pypospack.crystallography.Atoms`): a list 
            of atoms contained within the crystal structure.
        a0 (float): a scaling factor which scales the lattice vectors by a0.
            Default is 1.
        H (numpy.ndarray): a numpy array containing the lattice vectors as
            row vectors.  Default is [[1,0,0],[0,1,0],[0,0,1]]
        ptol (:obj:`float`,optional): tolerance in which to find an atom/layer.
            Defaults to None.
        vacancies (:obj:`list` of :obj:`pypospack.crystal.Atoms`): a
            list of vacancy sites contained within the crystal structure.
        interstitials (:obj:`list` of :obj:`pypospack.crystal.Atoms`): a
            list of interstitial sites contained within the crystal structure.

    """ 
    def __init__(self, obj=None):
        self.ptol = 1e-3

        if obj is None:
            self._noncopy_init()
        elif isinstance(obj,ase.atoms.Atoms):
            self._copy_init_ase(obj)
        elif isinstance(obj,SimulationCell):
            self._copy_init_pypospack(obj)
            
    def _noncopy_init(self):
        self.comment = "generated by pypospack"
        self.ptol = 1.e-5 #tolerance in which to find an atom
        self.a0  = 1.0
        self.H = np.array([[1,0,0],
                           [0,1,0],
                           [0,0,1]])
        self.atomic_basis = []
        self.vacancies = []
        self.interstitials = []

    def _copy_init_pypospack(self,obj):
        self.comment = obj.comment
        self.a0 = obj.a0
        self.H = np.array(obj.H)
        
        self.atomic_basis = copy.deepcopy(obj.atomic_basis)
        self.vacancies = copy.deepcopy(obj.vacancies)
        self.interstitials = copy.deepcopy(obj.interstitials)

    def _copy_init_ase(self,obj):
        self.comment = "generated by pypospack from ase"
        self.a0 = 1.0
        self.H = np.copy(obj.cell)
        self.atomic_basis = []
        for a in obj:
            symbol = a.symbol
            position = cartesian2direct(a.position,self.H)
            for i in range(3):
                if position[i] < 0:
                    position[i] = 1 + position[i]
            self.atomic_basis.append(
                    Atom(symbol=symbol,position=position))
        self.vacancies = []
        self.interstitials = []

    @property
    def h1(self):
        """numpy.ndarray: which is the a1 lattice vector"""
        return np.array(self.H[0,:])

    @h1.setter
    def h1(self, h1):
        self.H[0,:] = np.array(h1)

    @property
    def h2(self):
        """numpy.array: which is the a2 lattice vector"""
        return np.array(self.H[1,:])

    @h2.setter
    def h2(self,h2):
        self.H[1,:] = np.array(h2)
 
    @property
    def h3(self):
        """numpy.array: this is the a3 lattice vector"""
        return np.array(self.H[2,:])

    @h3.setter
    def h3(self,h3):
        self.H[2,:] = np.array(h3)

    @property
    def a1(self):
        """float: the length of h1 lattice vector"""
        a0 = self.a0
        h1 = self.h1
        a1 = a0 * h1.dot(h1)**0.5
        return a1

    @property
    def a2(self):
        """float: the length of the h2 lattice vector"""
        a0 = self.a0
        h2 = self.h2
        a2 = a0 * h2.dot(h2)**0.5
        return a2

    @property
    def a3(self):
        """float: the length of the h3 lattice vector"""
        a0 = self.a0
        h3 = self.h3
        a3 = a0*h3.dot(h3)**0.5
        return a3

    @property
    def b1(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)

        b1 = 2 * np.pi * np.cross(a2,a3) / np.dot(a1, np.cross(a2,a3))
        return b1

    @property
    def b2(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)


        b2 = 2 * np.pi * np.cross(a3,a1) / np.dot(a2, np.cross(a3,a1))
        return b2

    @property
    def b3(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = np.array(self.h1)
        a2 = np.array(self.h2)
        a3 = np.array(self.h3)

        b3 = 2 * np.pi * np.cross(a1,a2) / np.dot(a3, np.cross(a1,a2))
        return b3
        
    @property
    def n_atoms(self):
        """float: the number of atoms in the structure"""
        n_atoms = len(self.atomic_basis)
        return n_atoms

    @property
    def symbols(self):
        """list of str: a list of the symbols in the structure"""

        symbols = list(set([a.symbol for a in self.atomic_basis+self.interstitials]))
        return symbols

    def check_if_atom_exists_at_position(self,symbol,position):
        """determines if there is an atom at a position

        This code looks for atoms in the list of atoms in the atomic_basis.

        Returns:
            tuple: The first element if true in an atom exists at the location,
                False if it doesn't exist.   The second element is an int which
                indicates the index of the atom at that position.  The index is
                negative if it is the list of interstitials.
        """

        return_value = (False,None)
        # check to see if atom exists in the atomic basis
        for i,a in enumerate(self.atomic_basis):
            diff = [abs(position[i]-a.position[i]) for i in range(3)]
            if max(diff) < self.ptol:
                return_value = (True,i)

        return return_value

    def add_atom(self, symbol, position):
        """add an atom to the structure

        Checks to see if an existing atom exists at the position we are trying
        to add an atom, then if the position is empty.  The atom is added then
        added to the list of interstitials.

        Args:
            symbol (str): the symbol of the atom to be added
            position (str): the position of the atom to be added

        Raises:
            ValueError: If an atom already exists in the position.
        """
        (is_atom_exists, __) = self.check_if_atom_exists_at_position(symbol,position)
        if is_atom_exists:
            err_msg = "Tried to add {} @ {} an atom already there\n"
            err_msg = err_msg.format(symbol,position)
            err_msg += 'atomic basis:\n'
            for i,a in enumerate(self.atomic_basis):
                err_msg += ",".join([str(i),a,a.symbol,a.position])
            raise ValueError(err_msg)
        else:
            self.atomic_basis.append(Atom(symbol,position))

    def remove_atom(self,symbol, position):
        """ remove an atom from the structure

        This method checks for atom at the position, if an atom exists.  It is
        removed from the structure then the position is recorded as a vacancy.

        Args:
            symbol (str): the symbol of the atom
            position (:obj:`list` of :obj:`float`): the position of the atom

        """
        self.ptol = 1e-4
        for i,a in enumerate(self.atomic_basis):
            if (a.symbol == symbol):
                diff = [abs(position[j]-a.position[j]) for j in range(3)]
                is_atom = True
                for j in range(3):
                    if diff[j] >= self.ptol:
                        is_atom = False
                if is_atom:
                    self.atomic_basis.remove(self.atomic_basis[i])
                    return

        # no atom found
        err_msg = "Tried to remove {} @ {}, no atom found"
        err_msg = err_msg.format(symbol,position)
        raise ValueError(err_msg)

    def add_interstitial(self,symbol,position):
        """ add an interstitial to the atomic basis

        Args:
            symbol (str): the symbol of the atom
            position (:obj:`list` of :obj:`float`): the position of the atom
        """
        self.add_atom(symbol,position)
        self.interstitials.append([symbol,position])

    def add_vacancy(self,symbol,position):
        """ create a vacancy
        
        Creates a vacancy by removing an atom from the atomic basis

        Args:
            symbol (str): the symbol of the atom
            position (:obj:`list` of :obj:`float`): the position of the atom
        """
        self.remove_atom(symbol,position)
        self.vacancies.append([symbol,position])

    def get_number_of_atoms(self, symbol=None):
        if symbol is None:
            n_atoms = len(self._atoms)
        else:
            n_atoms = 0
            for atom in self.atomic_basis:
                if (atom.symbol == symbol):
                    n_atoms += 1
        return n_atoms

    def normalize_h_matrix(self):
        # change the h_matrix where the lattice parameter is 1.0
        for i in range(3):
            self.H[i,:] = self.H[i,:] * self.a0
        self.a0 = 1.0

        # calculate the norm of the h1 vector
        norm_h1 = self.a1
        # the norm of the h1 vector is the lattice parameter
        self.a0= norm_h1

        # normalize the h-matrix by the norm of h1
        for i in range(3):
            self.H[i,:] = self.H[i,:] / norm_h1

    def set_lattice_parameter(self, a0):
        # change the h_matrix where the lattice parameter is 1.0
        for i in range(3):
            self.H[i,:] = self.H[i,:] * self.a0

        self.a0 = a0
        for i in range(3):
            self.H[i,:] = self.H[i,:] / self.a0

    def __str__(self):
        str_out = "a = {}\n".format(self.a0)
        str_out += "atom list:\n"
        for a in self.atomic_basis:
            str_out = "{} {:10.6f} {:10.6f} {:10.6f} {:10.6f}"
            str_out += str_out.format(a.symbol,
                                      a.position[0],
                                      a.position[1],
                                      a.position[2],
                                      a.magnetic_moment)
 
def make_super_cell(structure, sc):
    """makes a supercell from a given cell

    Args:
        structure (pyflamestk.base.Structure): the base structure from which
            the supercell will be made from.
        sc (:obj:`list` of :obj:`int`): the number of repeat units in the h1, h2, and h3 
            directions
    """
    assert isinstance(structure,SimulationCell)

    supercell = SimulationCell()
    supercell.structure_comment = "{}x{}x{}".format(sc[0],sc[1],sc[2])

    # set lattice parameter
    supercell.a0 = structure.a0 

    # set h_matrix
    H = np.zeros(shape=[3,3])
    for i in range(3):
        H[i,:] = structure.H[i,:] * sc[i]
    supercell.H = H.copy()

    # add supercell atoms
    for i in range(sc[0]):
        for j in range(sc[1]):
            for k in range(sc[2]):
                for atom in structure.atomic_basis:
                    symbol = atom.symbol
                    position = atom.position
                    position = [(i+position[0])/sc[0],\
                                (j+position[1])/sc[1],\
                                (k+position[2])/sc[2]]
                    supercell.add_atom(symbol,position)

    # return a copy of the supercell
    return copy.deepcopy(supercell)

class RadialDistributionFunction(object):

    def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
        """Compute the three-dimensioanl pair correlation function
        
        Compute the three-dimensional pair correlation function for a set of
        spherical particles contained in a cube with side length S.  This simple
        function finds reference particles such that a sphere of radius rMax drawn
        around the particle will fit entirely within the cube, eliminating the need
        to compensate for edge effects.  If no such particles exist, an error is
        returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

        Arguments:
            x(numpy.ndarray): an array of x positions of centers of particles
            y(numpy.ndarray): an array of y positions of centers of particles
            z(numpy.ndarray): an array of z positions of centers of particles
            S(list):length of each side of the cube in space
            rMax(float): outer diameter of largest spherical shell
            dr(float): increment for increasing radius of spherical shell
        
        Returns:
            tuple: first element is a numpy array containing the correlation 
                function g(r) radii. the second element is a numpy array 
                containing the radii of the spherical shells used to compute 
                g(r).  the final element are the indices of the refence particles
        """
        from numpy import zeros, sqrt, where, pi, mean, arange, histogram

        # Find particles which are close enough to the cube center that a sphere of radius
        # rMax will not cross any face of the cube
        bools1 = x > rMax
        bools2 = x < (S - rMax)
        bools3 = y > rMax
        bools4 = y < (S - rMax)
        bools5 = z > rMax
        bools6 = z < (S - rMax)

        interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
        num_interior_particles = len(interior_indices)

        if num_interior_particles < 1:
            raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                    will lie entirely within a cube of side length S.  Decrease rMax\
                    or increase the size of the cube.")

        edges = arange(0., rMax + 1.1 * dr, dr)
        num_increments = len(edges) - 1
        g = zeros([num_interior_particles, num_increments])
        radii = zeros(num_increments)
        numberDensity = len(x) / S**3

        # Compute pairwise correlation for each interior particle
        for p in range(num_interior_particles):
            index = interior_indices[p]
            d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
            d[index] = 2 * rMax

            (result, bins) = histogram(d, bins=edges, normed=False)
            g[p,:] = result / numberDensity

        # Average g(r) for all interior particles and compute radii
        g_average = zeros(num_increments)
        for i in range(num_increments):
            radii[i] = (edges[i] + edges[i+1]) / 2.
            rOuter = edges[i + 1]
            rInner = edges[i]
            g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

        return (g_average, radii, interior_indices)
        # Number of particles in shell/total number of particles/volume of shell/number density
        # shell volume = 4/3*pi(r_outer**3-r_inner**3)

def get_fcc_nearest_neighbor_distance(a0,NN):
    assert isinstance(a0,float)
    assert isinstance(NN,int) or isinstance(NN,float)

    nn_distances = [
        0,
        0.707 * a0,
        1.000 * a0,
        1.225 * a0,
        1.414 * a0,
        1.581 * a0
    ]

    import math
    return nn_distances[math.ceil(NN)] \
            + (NN%1)*(nn_distances[math.floor(NN)]-nn_distances[math.ceil(NN)])
