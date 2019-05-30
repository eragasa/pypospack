__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2018"
__license__ = "Simplified BSD License"
__version__ = 20180524

import copy,inspect
import numpy as np
from collections import OrderedDict
from pypospack.potential import EamEmbeddingFunction

from pypospack.potential.eam_embedding_eos import \
        get_omega,\
        get_pair_energy_at_a,\
        get_density_at_a

def func_zopemishin_eos(a,a0,B0,E0,beta,lattice_type='fcc'):
    """ return the energy value from th equation of state based on the lattice parameter

    Args:
        a (float)(np.ndarray): an array of lattice parameter values to be evaluation
        a0 (float): the equilibrium lattice constant value, in Angstroms
        B0 (float): the bulk modulus value, GPA
        E0 (float): the cohesive energy at equilibrium, in eV
        beta (float): adjustable parameter, dimensionless
        lattice_type (str): the lattice name of the structure
    Returns:
        (float)(np.ndarray): the energy implied by the equation of state evaluated at the points
        implied by the argument a

    """
    x = a/a0-1 
    Omega0 = get_omega(a0)
    if False:
        print('a0:',a0)
        print('B0:',B0)
        print('E0:',E0)
        print('beta:',beta)
        print('lattice_type:',lattice_type)
        print('Omega0:',Omega0)
        print('quotient:',-(9*Omega0*B0)/E0)
    alpha = np.sqrt(-(9*Omega0*B0)/E0)
    E = E0 * (1+alpha*x+beta*(alpha**3)*(x**3)*(2*x+3)/((x+1)**2)) * np.exp(-alpha*x)
    return E

def func_zopemishin_embedding_function(rho,a0,B0,E0,beta,lattice_type='fcc'):
    """ fits the embedding function to the zope mishin equation of state

    This function only exists as a implementation prototype for optimization routines
    which requires an encapsulated function.  The EosFitter is a more general implementation,
    which should be used in most cases.

    Args:
        rho (numpy.ndarray): a list of densities for the embedding function to be evaluated at
        a0 (float): the equilibrium lattice parameters
        B0 (float): the equilibrium bulk modulus.
        beta (float): a shape parameter
        lattice_type (str): the type of lattice parameter, currently only 'fcc' is implemented

    Returns:
        (numpy.ndarray): evaluations of the embedding function for rho

    """

    func_density = o.formalisms['density']['Ni']
    func_density_param = o.parameters['popt']['density']['Ni']

    func_pair = o.formalisms['pair']['NiNi']
    func_pair_param = o.parameters['popt']['pair']['NiNi']

    func_eos = o.formalisms['embedding']['Ni']

    a_min=0
    a_max=10000
    a_tol=1.e-8

    if isinstance(rho,np.ndarray):
        rho_ = rho.tolist()
    else:
        rho_ = rho

    # calculate astars
    a = np.zeros(len(rho_))
    for i,rhostar in enumerate(rho_):
        a[i] = brentq(f=get_density_at_a,a=a_min,b=a_max,
                args=(rhostar,func_density,func_density_param),
                xtol=a_tol)
    
    E_eos = func_zopemishin_eos(a,a0,B0,E0,beta)
    E_pair = get_pair_energy_at_a(a,func_pair,func_pair_param)
    E_embedding = E_eos - E_pair

    return E_embedding

class ZopeMishinEosEmbeddingFunction(EamEmbeddingFunction):
    """Implementation of the Eam Embedding class for the universal equation of state
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
    
    def __init__(self,
            symbols,
            obj_density_function=None,
            obj_pair_function=None,
            lattice_type = None,
            lattice_a0 = None,
            parameters=None):

        EamEmbeddingEquationOfState.__init__(self,symbols=symbols)
        
        # define member variables, and initialize with None
        self.parameters =None
        self.obj_density_fn = None
        self.obj_pair_fn = None
        self.lattice_type = None
        self.lattice_a0 = None

        # process arguments
        if parameters is not None:
            self.parameters = copy.deepcopy(parameters)

        if obj_density_function is not None:
            self.obj_density_fn = obj_density_function
        
        if obj_pair_function is not None:
            self.obj_pair_fn = obj_pair_function

        if lattice_type is not None:
            self.lattice_type = lattice_type

        if lattice_a0 is not None:
            self.lattice_a0 = lattice_a0
    def __init__(self,
            symbols):
        self.embedding_func_parameters = ['F0','gamma','F1']
        EamEmbeddingFunction.__init__(self,
                symbols=symbols,
                potential_type = 'eam_embed_bjs')

    def _init_parameter_names(self):
        PARAMETER_NAME_FORMAT = "{s}_{pn}"
        PARAMETER_NAMES_EA_SYMBOL = ['ecoh','latticetype','B','a0']

        self.parameter_names = []

        for s in self.symbols:
            for pn in PARAMETER_NAMES_EA_SYMBOL:
                pn_string = PARAMETER_NAME_FORMAT.format(s=s,pn=pn)
                self.parameter_names.append(pn_string)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def get_embedding_energies(self, rho,
            func_eos, func_eos_param,
            func_density, func_density_param,
            func_pair, func_pair_param,
            a_min=0,a_max=10000,a_tol=1.e-8):

            rho_ = rho.tolist()

            a = np.zeros(len(rho_)) # pre-allocate
            for i,rhostar in enumerate(rho_):
                    a[i] = brentq(f=get_density_at_a,
                        args=(rhostar,func_density,func,func_density_param),
                        a=a_min,b=a_max,xtol=a_tol)
            arg_names = [
                k for k in  inspect.getargspec(func_eos)[0] \
                    if k not in ['rho,''lattice_type']
            ]
            args = [func_eos_param[k] for k in arg_names]
            E_eos = func_eos(a,*args) # energy from equation of state
            E_pair = get_pair_energy_at_a(a,func_pair,func_pair_param)

            # assume that E_eos = E_embedding + E_pair, therefore
            E_embedding = E_eos - E_pair

    def evaluate(self,
            rho,
            r,
            parameters,
            o_pair=None,
            o_density=None,
            a_min=0.,
            a_max=10000.,
            a_tol=1.e-8):

        # need to add error handling code here
        if o_pair is not None:
            self.obj_pair_fn = o_pair
        else:
            if self.obj_pair_fn is None:
                raise ValueError("pair potential function attribute is None")
        
        # need to add error handling code here as well
        if o_density is not None:
            self.obj_density_fn = o_density
        else:
            if self.obj_density_fn is None:
                raise ValueError("electron density function is None")

        assert isinstance(parameters,OrderedDict)
        # we can define rho at a point as a function
       
        # nested function for development
        def rhofxn(a,args):
            """
            a - the lattice parameter for which we are looking for the total electron density
            args[0] - parameters
            args[1] = rho bar.
            """
            rho_fxn_parameters = OrderedDict()
            _rho_bar= args[1]
            _s = s

            for k,v in args[0].items():
                if k.startswith('d_{}'.format(_s)):
                    rho_fxn_parameters[k[2:]] = v

            lattice_type = args[0]['e_{}_latticetype'.format(_s)]
            if lattice_type == 'fcc':
                n_NN = [12,6,24,12,24,8]
                d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
            
            _rho = 0.

            if isinstance(r,np.ndarray) and isinstance(rho,np.ndarray):
                _density = self.obj_density_fn.density_evaluations[s]
            
                for i in range(len(n_NN)):
                    _rho += n_NN[i] * np.interp(
                            d_NN[i],
                            r,
                            self.obj_density_fn.density_evaluations[s]
                        )

            else:
                for i in range(len(n_NN)):
                    if len(rho_fxn_parameters) == 0:
                        _rho += n_NN[i] * self.obj_density_fn.evaluate(
                                d_NN[i],
                                rho_fxn_parameters
                            )[s]
                    else:
                        _rho += n_NN[i] * self.obj_density_fn.evaluate(
                                d_NN[i],
                                rho_fxn_parameters
                            )[s]
            
            return _rho-_rho_bar
        
        def pairfxn(r,args,pair_key=None,r_vector=None):

            pair_fxn_parameters = OrderedDict()
            for k,v in args.items():
                if k.startswith('p_{}{}'.format(s,s)):
                    pair_fxn_parameters[k[2:]] = v

            lattice_type = args['e_{}_latticetype'.format(s)]
            if lattice_type == 'fcc':
                n_NN = [12,6,24,12,24,8]
                d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
            
            _pot = 0.

            if isinstance(r_vector,np.ndarray):
                _pair_key = '{}{}'.format(s,s)
                for i in range(len(n_NN)):
                    _pot += n_NN[i] * np.interp(
                            d_NN[i],
                            r_vector,
                            self.obj_pair_fn.potential_evaluations[_pair_key]
                        )
            else:
                print('notinterp')
                for i in range(len(n_NN)):
                    if len(pair_fxn_parameters) == 0:
                        _pot += n_NN[i] * self.obj_pair_fn.evaluate(
                                d_NN[i],
                                pair_fxn_parameters
                            )['{}{}'.format(s,s)]
                    else:
                        _pot += n_NN[i] * self.obj_pair_fn.evaluate(
                                d_NN[i],
                                pair_fxn_parameters
                            )['{}{}'.format(s,s)]

            # ?? remove double counting of pairs
            _pot = 0.5 * _pot
            return _pot

        self.embedding_evaluations = OrderedDict()
        #unpack variable array.
        for s in self.symbols:
            embed_vals = np.empty_like(rho)
            try:
                _ecoh = parameters['e_{}_ecoh'.format(s)]
            except KeyError as e:
                _ecoh = parameters['{}_ecoh'.format(s)]

            try:
                _latticetype = parameters['e_{}_latticetype'.format(s)]
            except KeyError as e:
                _latticetype = parameters['{}_latticetype'.format(s)]

            try:
                _B = parameters['e_{}_B'.format(s)]
            except KeyError as e:
                _B = parameters['{}_B'.format(s)]
            
            try:
                _a0 = parameters['e_{}_a0'.format(s)]
            except KeyError as e:
                _a0 = parameters['{}_a0'.format(s)]

            if _latticetype == 'fcc':
                _r0 = _a0/np.sqrt(2.)
                _n_atoms_per_unit_cell = 4

            for k,rhostar in enumerate(rho):
        
                # for each rho, we need to determine the lattice parameter which produces
                # the rho for an atom in the corresponding lattice type.
                try:
                    a = brentq(
                            f=rhofxn,
                            a=a_min,
                            b=a_max,
                            args=[parameters,rhostar,s,r],
                            xtol=a_tol)
                except ValueError as e:
                    # this error is thrown due to the brent bounding bracket [a, b] not including a zero
                    raise PypospackBadEamEosError(parameters=parameters)

                # here we determine astar as in equation 5 in Foiles. Phys Rev B (33) 12. Jun 1986 
                # 160.22 is the conversion with _esub in eV, _B in GPa, and _omega in Angs^3
                _esub = abs(_ecoh) 
                _omega = (_a0**3)/_n_atoms_per_unit_cell 
                astar = ((a/_a0)-1)/((_esub/(9*_B*_omega) * 160.22)**0.5)

                # now determine the rose energy as in Equation 4
                _e_rose = -_esub*(1+astar)*np.exp(-astar)
                
                # now find the pair potential for the lattice constant
                _pair_key = '{}{}'.format(s,s)
                _e_pot = pairfxn(a,parameters,_pair_key,r)

                # now the _e_rose is the total energy of the system
                # e_rose = e_pot + e_embedding
                # e_embedding = e_rose - e_pot
                embed_vals[k] = _e_rose - _e_pot

            self.embedding_evaluations[s] = copy.deepcopy(embed_vals)

    def evaluate(self,rho,parameters):
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
            rho_cut(float): This would be the density cutoff.  However the
                embedding energy is increasing with increasing electron
                density so the a r_cut has no physical meaning.  Any 
                variable passed into r_cut will be ignored.
        """
        # attribute.parameters[p] <--- arg:parameters[p]
        for s in self.symbols:
            for p in self.embedding_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameters[pn] = parameters[pn]

        self.embedding_evaluations = OrderedDict()
        for s in self.symbols:

            # get parameters
            F0 = self.parameters['{}_F0'.format(s)]
            gamma = self.parameters['{}_gamma'.format(s)]
            F1 = self.parameters['{}_F1'.format(s)]
            
            with np.errstate(all='raise'):
                self.embedding_evaluations[s] \
                    = F0*(1-gamma*np.log(rho))*rho**gamma + F1*gamma
        
        return copy.deepcopy(self.embedding_evaluations)
