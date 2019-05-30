import numpy as np
from pypospack.potential.eamdens_mishin2003 import func_cutoff_mishin2003
from pypospack.potential.eamdens_mishin2003 import func_mishin2003_density
from pypopsack.potential.eamdens_mishin2003 import func_mishin2003_density_w_cutoff

def create_r(rmax,n):
    r = r_max/n * np.linspace(1,n,n)
    return r

if __name__ == "__main__":
    from collections import OrderedDict
    import matplotlib.pyplot as plt

    r_max = 3.
    r_N = 1000
    fig, ax = plt.subplots(1,1)

    density_mishin2003_parameters = OrderedDict([
            ('r0',1.),
            ('A0',1.),
            ('B0',1.),
            ('C0',1.),
            ('y',1.),
            ('gamma',1.),
            ('rc',2.),
            ('h',2.)
        ])
    r = create_r(rmax=r_max,n=r_N)
    dens_w_cutoff = func_mishin2003_density_w_cutoff(r=r,**density_mishin2003_parameters)
    ax.plot(r,dens_w_cutoff)
    plt.show()

