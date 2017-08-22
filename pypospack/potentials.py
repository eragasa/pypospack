import numpy as np
import scipy.optimize

class PyFlamesTkError(Exception):
    def __init__(self,value):
        self.value = value

    def __str__(self):
        return self.value
        
class Potential(object):
    """
    Class: Potential
    Author: Eugene J. Ragasa
    Date: 2/26/2016
    
    This is class is an abstract class for empirical potentials.

    Potential.__init__ must    
    """
    def __init__(self, name, pot_name, atoms):
        self.name = name
        self.potential_name = pot_name
        self.atoms = {}
        for atom in atoms:
            self.atoms[atom] = {}
        self.parameter_list = {}
        
    def evaluate(self):
        pass
    def set_parameter_list(self):
        pass
    def get_parameter_list(self):
        pass

class Buckingham(Potential):
    def __init__(self,name,atoms):
        Potential.__init__(self, name, 'buckingham', atoms)
        for atom in self.atoms:
            print('atom',atom)
            
        
class EamPotential(Potential):
    def __init__(self, name,cutoff, embed_type, dens_type, pair_type, atoms):
        Potential.__init__(name, 'eam/analytical',atoms)
        self.cutoff = 0
        self.embed_type = embed_type
        self.dens_type  = dens_type
        self.pair_type  = pair_type
        self._initialize_potentials()
        
    def _initialize_potentials():
        self.embedding = {}
        self.density   = {}
        self.pair      = {}

        for atom_i in self.atoms:
            self.embedding[atom_i] = []
            self.density[atom_i]   = []

        for atom_i in self.atoms:
            for atom_j in self.atoms:
                pair_name     = "{}{}".format(atom_i,atom_j)
                inv_pair_name = "{}{}".format(atom_j,atom_i)
                if not(inv_pair_name in self.atoms):
                    self.pair[pair_name] = []
            
def func_psi(r):
    ''' Implements cutoff function for smoothly bringing function to zero '''
    if type(r) == np.ndarray:
        s = np.empty_like(r)
        for i in range(len(r)): 
            x = r[i]
            if x > 1.0:
                s[i]= 0.
            elif ((x > 0.0) and (x <= 1.0)):
                s[i] = ( -6.*x**5 + 15.*x**4 - 10.*x**3 + 1.)
            else:
                s[i] = 1.
    else:
        if r > 1.0:
            s = 0.
        elif ((r > 0.0) and (r <= 1.0)):
            s = ( -6.*r**5 + 15.*r**4 - 10.*r**3 + 1.)
        else:
            s = 1.

    return s

def func_exponential_density(r,a):
    val = 0
    return val

def func_morse(r,De,am,r0,rp, cutoff = 0):
    pairvals = []
    if cutoff != 0:
        pairvals = De * ( ( 1.- np.exp( -am*(r-r0) ) )**2. - 1. ) * func_psi( (r-rp)/(cutoff-rp) )
    else:
        pairvals = De * ( ( 1.- np.exp( -am*(r-r0) ) )**2. - 1. )
    return pairvals
    
def pairpotential_voterchen(r,dm,rm,am,rcut):
    # Ref: Chen, Srolovitz, Voter, J.Mater.Res. 4,1, 1989
    '''
    r  - radius
    dm - (Dm) the depth
    rm - (r_m) distance to the minimum
    am - (alpha_m) measure of curvature at min of Morse pot
    '''
    fcut = np.exp(1/(r-rcut))
    phi  = (dm *(1-np.exp(-am*(r-rm)))**2-dm)*fcut
    return phi
    
def electrondensity_voterchen(r,c1,c2,rcut):
    # Ref: Chen, Srolovitz, Voter, J.Mater.Res. 4,1, 1989
    fcut = np.exp(1/(r-rcut))
    
    rho = c1*(r**6)*(np.exp(-c2*r)+2**9*np.exp*(-2*c2*r))*fcut
    return rho

def cuttoff_exponential(r,r_cut):
    fcut = np.exp(1/(r-rcut))
    return fcut

def rhofxn(a,rho0,r0,lambda0,rd,rhostar):
    ''' calculates ideal e- density based on exponential functional form 
        data input format:  rho0  r0    lambda0  rd   rhostar ''' 
    return rho0*(12.*exp(-(a/sqrt(2.)-r0)/lambda0)  * psi( (a/sqrt(2.)-rd) / (globalcutoff-rd) )
               +  6.*exp(-(a-r0)/lambda0)           * psi( (a-rd)          / (globalcutoff-rd) )
               + 24.*exp(-(a*sqrt(1.5)-r0)/lambda0) * psi( (a*sqrt(1.5)-rd) / (globalcutoff-rd) )
               + 12.*exp(-(a*sqrt(2.)-r0)/lambda0)  * psi( (a*sqrt(2.)-rd) / (globalcutoff-rd) )
               + 24.*exp(-(a*sqrt(2.5)-r0)/lambda0) * psi( (a*sqrt(2.5)-rd) / (globalcutoff-rd) )
               +  8.*exp(-(a*sqrt(3.)-r0)/lambda0)  * psi( (a*sqrt(3.)-rd) / (globalcutoff-rd) )
            ) - rhostar

def embedding_eos_voterchen(rho,ecoh,am,r0,rho0,lambda0,lambdarose,de,rp,rd):
    ''' 
    implements Foiles-style embeding function where the density function and
    pair potential marches the Rose equation of state, for the voter-srolovitz-voter
    electron density
    parameter:
    rho  -
    ecoh - cohesive energy or sublimination energy
    a    - equilibrium lattice parameter
    r0   - 
    rho0 - 
    lambda0
    lambdarose
    de
    rp
    rd
    '''
    embedvals = np.empty_like(rho)
    k=0
    for rho_i in rho:
        #solve density for lattice constant (a) where density (rhostar) is found
        
        # Uses Brent(1973) method for finding the zero of a function
        rho_p = (rho0,r0,lambda0,rd,rhp_i)        
        a = scipy.optimize.brentq(electrondensity_voterchen,
                                  0.,
                                  10000.,
                                  rho_p,
                                  1.0e-8) #lattice constant where rhostar is found
        #find E_Rose for lattice constants
        a_star = (a/a0 - 1.)/(E_coh/9 * bulk_mod * atm_vol)
        astar = (a-r0*np.sqrt(2.)) / (lambdarose*np.sqrt(2.)*r0) 
        Erose = Ecoh*(1+astar)*exp(-astar) 
        #find pair potential for lattice constant
        pp = (De,am,r0,rp)
        Epot = 12.*fpair(a/sqrt(2.),pp)+6.*fpair(a,pp)+24.*fpair(sqrt(1.5)*a,pp)+12.*fpair(sqrt(2.)*a,pp)+24.*fpair(sqrt(2.5)*a,pp)+8.*fpair(sqrt(3.)*a,pp)
        #calculate value of embedding fxn
        embedvals[k] = Erose - 0.5*Epot
        k += 1

    return embedvals
    
