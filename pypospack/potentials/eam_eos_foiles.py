import numpy as np
from scipy.optimize import brentq


class EamEmbeddingEquationOfState(object):

    def __init__(self,parameters):
        self.density_fn = None
        self.pair_fn = None

    @property
    def density_function(self):
        return self.density_fn

    @property
    def pair_function(self):
        return self.pair_fn


    def equation_of_state(self,rho,parameters):
        if latt_type is 'fcc':

    def evaluate(rho,parameters,o_pair,o_density):
        pass

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
        F_xtol = 1.0e-8,
        F_xtol=1.0e-8):

    F_evals = np.empty_like(rho)


    for rhostar in rho:
        p_embedding = (rho0,r0,lambda0,rd,rhostar)
        astar = brentq(
                rhofxn,
                a=F_min,
                b=F_max,
                p_embedding,
                xtol=F_xtol
        )
