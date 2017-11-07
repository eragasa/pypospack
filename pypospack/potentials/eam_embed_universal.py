__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
from collections import OrderedDict
from pypospack.potential import EamEmbeddingFunction

class UniversalEmbeddingFunction(EamEmbeddingFunction):
    """
    Args:
        symbols(list of str)
    Attributes:
        symbols(list of str)
        potential_type(str): This is set to 'eamembed_universal'
        parameter_names(list of str)
        parameters(OrderedDict): The key is the symbol associated with the
            embedding function.  On initialization, the value of each parameter
            is set to None.
        embedding(OrderedDict): The key is the symbol associated with the 
            embedding function.
        N_rho(int)
        d_rho(float)
        rho_max(float)
        rho(numpy.ndarray)
    """    
    def __init__(self,symbols):
        self.embedding_func_parameters = ['F0','p','q','F1']
        EamEmbeddingFunction.__init__(self,
                symbols=symbols,
                potential_type='eam_embed_univeral')

    def _init_parameter_names(self):
        self.parameter_names = []
        for s in self.symbols:
            for p in self.embedding_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameter_names.append(pn)
    
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def evaluate(self,rho,parameters,r_cut=None):
        """

        Given a vector of electron densities, rho, passed in as variable
        r, and the associated parameters of the potential.  This method
        sets the embedding attribute.

        Args:
            rho(numpy.ndarray): This should be named as rho because it
                represents the electron density being evaluated.
            parameters(OrderedDict): This is a dictionary of the parameters
                of the embedding function for each atom.  The key is a
                string containing the ISO chemical symbol of the element.
                The value should be a numeric value.
            r_cut(float): This would be the density cutoff.  However the
                embedding energy is increasing with increasing electron
                density so the a r_cut has no physical meaning.  Any 
                variable passed into r_cut will be ignored.
        """

        # attribute.parameters[p] <--- arg:parameters[p]
        for s in self.symbols:
            for p in self.embedding_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameters[pn] = parameters[p]

        self.embedding = OrderedDict()
        for s in self.symbols:
            F0 = self.parameters['{}_F0'.format(s)]
            q = self.parameters['{}_q'.format(s)]
            p = self.parameters['{}_p'.format(s)]
            F1 = self.parameters['{}_F1'.format(s)]
            self.embedding[s] \
                    = F0*(q/(q-p)*rho**p-p/(q-p)*rho**q)+F1*rho
        return copy.deepcopy(self.embedding)

