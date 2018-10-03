import copy
from collections import OrderedDict
import numpy as np
from pypospack.potential import EamEmbeddingEquationOfState
from scipy.optimize import brentq


"""
This class implements the determination of the embedding function by inverse solution of the Rose-Vinet Equation of State[1].
This method is outlined in by Foiles.

References:
[1] Foiles, Baskes, and Daw., Phys. Rev. B., 1986
[2] Brent, R. P., Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4
"""

#def rhofxn(a,rho0,r0,lambda0,rd,rhostar):
#    ''' calculates ideal e- density based on exponential functional form 
#        data input format:  rho0  r0    lambda0  rd   rhostar ''' 
#    return rho0*(
#              12.*exp(-(a/sqrt(2.)-r0)/lambda0)  * psi( (a/sqrt(2.)-rd) / (globalcutoff-rd) )
#            + 6. *exp(-(a-r0)/lambda0)           * psi( (a-rd) / (globalcutoff-rd) )
#            + 24.*exp(-(a*sqrt(1.5)-r0)/lambda0) * psi( (a*sqrt(1.5)-rd) / (globalcutoff-rd) )
#            + 12.*exp(-(a*sqrt(2.)-r0)/lambda0)  * psi( (a*sqrt(2.)-rd) / (globalcutoff-rd) )
#            + 24.*exp(-(a*sqrt(2.5)-r0)/lambda0) * psi( (a*sqrt(2.5)-rd) / (globalcutoff-rd) )
#            + 8. *exp(-(a*sqrt(3.)-r0)/lambda0)  * psi( (a*sqrt(3.)-rd) / (globalcutoff-rd) )
#            ) - rhostar

def rose_equation_of_state(a,a0,e_coh):
    """
    Args:
        B0 - isothermal bulk modulus
        dB0_dP - first derivative of pressure w.r.t to Pressure
        omega = 
    """
    
    a_star = (a/a0-1)/sqrt
    return E

class RoseEquationOfStateEmbeddingFunction(EamEmbeddingEquationOfState):

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

    def rose_equation_of_state(self,s):
        cubic_lattice_types = "fcc,bcc,dia"

        # the equation of state function
        e_coh = self.parameters["{}_ecoh".format(s)]
        bulk_modulus = self.parameters["{}_B".format(s)]
        
        a0 = self.parameters["{}_a0"]
       
        # this should be replaced with ore general code which implements a base class
        # lattice to retrieve all this reference data

        if self.parameters["{}_latticetype"] is 'fcc':
            V = a0**3
            n_atoms_per_unit_cell = 4
            omega = V/n_atoms_per_unit_cell
        
        
        e_rose = - e_coh
        

    def equation_of_state(self,rho,parameters=None):
        if parameters is None:
            _p = self.parameters
        else:
            _p = parameters

        if latt_type is 'fcc':
            pass

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
                a = brentq(
                        f=rhofxn,
                        a=a_min,
                        b=a_max,
                        args=[parameters,rhostar,s,r],
                        xtol=a_tol)

                # here we determine astar as in equation 5 in Foiles. Phys Rev B (33) 12. Jun 1986 
                _esub = abs(_ecoh) * 1.602e-19   # eV to J conversions as well
                _omega = (_a0**3)/_n_atoms_per_unit_cell * (1e-10)**3
                astar = ((a/_a0)-1)/((_esub/(9*_B*_omega))**0.5)

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
if __name__ == "__main__":
    from collections import OrderedDict
    p = OrderedDict()

    # testing the constructor
    o = EamEmbeddingEquationOfState(parameters=p)
    assert o.density_fn is None
    assert o.pair_fn is None
    assert o.r_cut is None

    
