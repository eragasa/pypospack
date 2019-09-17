# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
from pypospack.exceptions import BadParameterException
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018,2019"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from collections import OrderedDict

class Potential(object):
    """base class for potential

    Args:
        symbols(list of str): a list of symbols to use for the potential
        potential_type(str): a description for the type of potential this is

    """

    # pylint: disable=too-many-instance-attributes
    PYPOSPACK_CHRG_FORMAT = "chrg_{s}"
    PYPOSPACK_PAIR_FORMAT = "{s1}{s2}_{p}"
    PYPOSPACK_3BODY_FORMAT = "{s1}{s2}{s3}_{p}"
    def __init__(self,
                 symbols,
                 potential_type=None,
                 is_charge=None):

        # define formatting strings

        self.potential = None
        self.symbols = symbols
        self.potential_type = potential_type
        self.is_charge = is_charge

        # these attributes will be initialized by _init_parameter_names
        self.symbol_pairs = None
        self.parameter_names = None
        self._init_parameter_names()

        # these attributes will be initialized by _init_parameter_names
        self.parameters = None
        self._init_parameters()

        # deprecated parameters here
        self.param = {}
        self.param_names = None         # list of str

    def _init_parameter_names(self):
        raise NotImplementedError

    def _init_parameters(self):
        raise NotImplementedError

    def evaluate(self, r, parameters, r_cut=False):
        """evaluate the potential

        Args:
            r(numpy.ndarray): a numpy array of interatomic distances
            parameters(OrderedDict): an dictionary of parameter values and keys
            r_cut(float,optional): the global cutoff for the potential
        """
        raise NotImplementedError

    def write_lammps_potential_file(self):
        """writes the lammps_potential file

        This method exists to write the lammps potential file to disk.  This 
        method needs to be overriden to be implemented, in classes that
        inherit from this class.
        """

        raise NotImplementedError

    def lammps_potential_section_to_string(self):
        """generates the lammps string for the lammps potential sections

        Returns:
            str: the string value for the LAMMPS potential section that goes into the potential.mod file.
        """

        raise NotImplementedError


    def write_gulp_potential_section(self):
        """writes the gulp potential for a GULP string"""
        raise NotImplementedError

    def gulp_potential_section_to_string(self):
        """generates the potential section for a GULP string"""
        raise NotImplementedError

    def _get_mass(self, element):
        amu = OrderedDict([
            ('B',  10.811),
            ('C',  12.0107),
            ('N',  14.0067),
            ('O',  15.999),
            ('Mg', 24.305),
            ('Al', 26.982),
            ('Si', 28.0855),
            ('Ni', 58.6934)
            ])

        try:
            return amu[element]
        except NameError as e:
            raise

    def _get_name(self, element):
        element_name = OrderedDict([
            ('B',  'boron'),
            ('C',  'carbon'),
            ('N',  'nitrogen'),
            ('O',  'oxygen'),
            ('Mg', 'magnesium'),
            ('Al', 'aluminum'),
            ('Si', 'silicon'),
            ('Ni', 'nickel')
            ])

        try:
            return element_name['B']
        except NameError as e:
            raise

if __name__ == "__main__":
    o = Potential(symbols=['Ni'])

