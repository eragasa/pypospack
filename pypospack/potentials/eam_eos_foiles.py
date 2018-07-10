import numpy as np
from pypospack.potential import EamEmbeddingEquationOfState
from scipy.optimize import brentq

class RoseEquationOfStateEmbeddingFunction(EamEmbeddingEquationOfState):

    def __init__(self,parameters):
        EamEmbeddingEquationOfState.__init__(self,parameters)

    def rose_equation_of_state(self):
        E = None
        return E

    def equation_of_state(self,rho,parameters=None):
        if parameters is None:
            _p = self.parameters
        else:
            _p = parameters

        if latt_type is 'fcc':
            pass

    def evaluate(rho,parameters,o_pair,o_density):
        embed_vals = np.empty_like(rho)
        _ecoh = parameters['ecoh']
        

def fembedFoiles(rho,params):
    ''' implements Foiles-style embeding function 
        (i.e. density and pair potential forced to match Rose EOS)
        parameter list:
        p[0]   p[1]     p[2] p[3] p[4]    p[5]       p[6]   p[7] p[8] 
        E_coh  a(morse) r0  rho0 lambda0 lambdarose De      rp   rd  '''
    Ecoh,am,r0,rho0,lambda0,lambdarose,De,rp,rd = params
    embedvals = empty_like(rho)
    k=0
    for rhostar in rho:
        #solve density for lattice constant (a) where density (rhostar) is found
        rhop = (rho0,r0,lambda0,rd,rhostar)
        a = brentq(rhofxn,0.,10000.,rhop,xtol=1.0e-8) #lattice constant where rhostar is found
        #find E_Rose for lattice constant
        astar = (a-r0*sqrt(2.)) / (lambdarose*sqrt(2.)*r0) 
        Erose = Ecoh*(1+astar)*exp(-astar) 
        #find pair potential for lattice constant
        pp = (De,am,r0,rp)
        Epot = 12.*fpair(a/sqrt(2.),pp)+6.*fpair(a,pp)+24.*fpair(sqrt(1.5)*a,pp)+12.*fpair(sqrt(2.)*a,pp)+24.*fpair(sqrt(2.5)*a,pp)+8.*fpair(sqrt(3.)*a,pp)
        #calculate value of embedding fxn
        embedvals[k] = Erose - 0.5*Epot
        k += 1

    return embedvals

def psi(r):
    ''' Implements cutoff function for smoothly bringing function to zero '''
    if type(r) == ndarray:
        s = np.empty_like(r)
        for i in xrange(len(r)): 
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

# This is the exponential density function at rho
def rhofxn(a,rho0,r0,lambda0,rd,rhostar):
    ''' calculates ideal e- density based on exponential functional form 
        data input format:  rho0  r0    lambda0  rd   rhostar ''' 
    return rho0*(
              12.*exp(-(a/sqrt(2.)-r0)/lambda0)  * psi( (a/sqrt(2.)-rd) / (globalcutoff-rd) )
            + 6. *exp(-(a-r0)/lambda0)           * psi( (a-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(1.5)-r0)/lambda0) * psi( (a*sqrt(1.5)-rd) / (globalcutoff-rd) )
            + 12.*exp(-(a*sqrt(2.)-r0)/lambda0)  * psi( (a*sqrt(2.)-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(2.5)-r0)/lambda0) * psi( (a*sqrt(2.5)-rd) / (globalcutoff-rd) )
            + 8. *exp(-(a*sqrt(3.)-r0)/lambda0)  * psi( (a*sqrt(3.)-rd) / (globalcutoff-rd) )
            ) - rhostar

def func_eam_embed_foiles(
        rho,
        E0,
        am,
        r0,
        rho0,
        lambda0,
        lambdarose,
        De,
        rp,
        rd,
        F_min = 0,
        F_max = 10000,
        F_xtol=1.0e-8):

    F_evals = np.empty_like(rho)


    for rhostar in rho:
        p_embedding = (rho0,r0,lambda0,rd,rhostar)
        astar = brentq(
                rhofxn,
                a=F_min,
                b=F_max,
                args=p_embedding,
                xtol=F_xtol
        )

if __name__ == "__main__":
    from pypospack.potential import 
    from collections import OrderedDict
    p = OrderedDict()

    # testing the constructor
    o = EamEmbeddingEquationOfState(parameters=p)
    assert o.density_fn is None
    assert o.pair_fn is None
    assert o.r_cut is None

    
