import inspect
from collections import OrderedDict

from pypospack.potential import EamEmbeddingFunction

def get_omega(a0,lattice_type='fcc'):
    V=a0**3
    if lattice_type == 'fcc':
        n_atoms = 4
    else:
        m = '{} is an unsupported lattice type'
        raise ValueError(m)

    return V/n_atoms

def get_density_at_a(a,
        rho_bar,
        func_density,
        func_density_param,
        lattice_type='fcc'):

    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')
        raise ValueError(m)

    #calculate total density
    if isinstance(a,np.ndarray):
        rho = np.zeros(len(r))
    else:
        rho = 0.
    arg_names = [k for k in inspect.getargspec(func_density)[0] if k != 'r']
    args = [func_density_param[k] for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        rho += n * func_density(r,*args)

    return rho - rhobar

def get_pair_energy_at_a(a,
        func_pair,
        func_pair_param,
        lattice_type='fcc'):


    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')

    # calculate total energy
    if isinstance(a,np.ndarray):
        phi = np.zeros(len(a))
    else:
        phi = 0.
    arg_names =[k for k in inspect.getargspec(func_pair)[0] if k != 'r']
    args = [func_pair_param[k] for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        phi += n * func_pair(r,*args)

    # divide by 2, to eliminate double counting
    return 0.5 * phi


class EamEmbeddingEquationOfState(EamEmbeddingFunction):

    potential_type = 'eam_eos_base'
    def __init__(self,
                 symbols,
                 obj_density_function=None,
                 obj_pair_function=None,
                 lattice_type=None,
                 lattice_a0=None,
                 parameters=None
                 ):

        potential_type_ = 'eam_eos_base'

        EamEmbeddingFunction.__init__(self,
                                      symbols=symbols,
                                      potential_type='eam_eos_base')

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
