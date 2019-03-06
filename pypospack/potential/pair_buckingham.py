__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyrirght (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

class BuckinghamPotential(PairPotential):
    """ Implementation of the Buckingham Potential

    This class provides an interface for the management of parameters
    to and from different molecular dynamics and lattice dynamics programs.

    Args:
        symbols(list): a list of chemicals
    """

    def __init__(self,symbols):
        self.pair_potential_parameters = ['A','rho','C']
        PairPotential.__init__(self,
                symbols,
                potential_type='buckingham',
                is_charge=True)

    def _init_parameter_names(self):
        self.symbol_pairs = list(determine_symbol_pairs(self.symbols))
        self.parameter_names = []

        for s in self.symbols:
            self.parameter_names.append(
                self.PYPOSPACK_CHRG_FORMAT.format(s=s))

        for sp in self.symbol_pairs:
            for p in self.pair_potential_parameters:
                self.parameter_names.append(
                        self.PYPOSPACK_PAIR_FORMAT.format(
                            s1=sp[0],s2=sp[1],p=p))

        return list(self.parameter_names)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    def evaluate(self,r,parameters,r_cut=False):
        raise NotImplementedError

    # same as parent class
    def lammps_potential_section_to_string(self,parameters=None,rcut=10.0):
        """get the string for the lammps potential section

        Args:
            parameters(OrderedDict): a dictionary of parameter name keys with the associated values
            rcut(float): the global cutoff
        Returns:
            str: the string of the LAMMPS potential section
        """

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        # set masses
        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        # set groups
        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        # set chrg
        for i,s in enumerate(self.symbols):
            charge = self.parameters['chrg_{}'.format(s)]
            str_out += "set group {} charge {}\n".format(s,charge)
        str_out += "\n"

        str_out += 'variable R_cut equal {}\n'.format(rcut)
        str_out += '\n'
        str_out += 'pair_style buck/coul/long ${R_cut}\n'

        # set param values
        for i,si in enumerate(self.symbols):
            for j,sj in enumerate(self.symbols):
                if i <= j:
                    try:
                        A = self.parameters['{}{}_A'.format(si,sj)]
                        rho = self.parameters['{}{}_rho'.format(si,sj)]
                        C = self.parameters['{}{}_C'.format(si,sj)]
                        str_out += "pair_coeff {} {} {} {} {} {}\n".format(\
                                i+1,j+1,A,rho,C,'${R_cut}')
                    except KeyError as ke:
                        s = str(ke)
                        print('key_requested:',s)
                        print('keys:',self.parameters.keys())
                        raise
                    except TypeError as te: 
                        s = str(te)
                        print(self.param_dict)
                        print('key_requested:',s)
                        print('keys:',self.parameters.keys())
                        raise
        

        return str_out
    
    # overrides the parents class
    def gulp_potential_section_to_string(self,parameters=None,r_cut=10.0):
        """ get GULP potential to string

        The buckingham potential is a charged potential and so the charges
        associated with the potential are also part of the potential."
        """
        if parameters is not None:
            for pn in self.parameters:
                self.parameters[pn] = parameters[pn]

        str_out = 'species\n'
        for s in self.symbols:
            chrg=self.parameters['chrg_{}'.format(s)]
            str_out += "{s} core {chrg}\n".format(s=s,chrg=chrg)

        str_out += 'buck\n'

        for symbols in self.symbol_pairs:
            s1 = symbols[0]
            s2 = symbols[1]
            sp = "{}{}".format(s1,s2)

            # get parameters
            A = parameters['{}_A'.format(sp)]
            rho = parameters['{}_rho'.format(sp)]
            C = parameters['{}_C'.format(sp)]

            str_out += "{s1} core {s2} core {A} {rho} {C} {r_cut}\n".format(
                    s1=s1,s2=s2,A=A,rho=rho,C=C,r_cut=r_cut)
        
        return str_out
    
    # same as parent class
    def phonts_potential_section_to_string(self):
        raise NotImplementedError 
    
    # same as parent class
    def write_lammps_potential_file(self):
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        raise NotImplementedError
