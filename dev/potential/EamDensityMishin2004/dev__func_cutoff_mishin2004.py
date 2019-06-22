import numpy as np
from eamdens_mishin2004 import func_cutoff_mishin2004

def create_r(rmax,n):
    r = r_max/n * np.linspace(1,n,n)
    return r

if __name__ == "__main__":
    from collections import OrderedDict
   
    r_max = 3.
    r_N = 1000

    cutoff_mishin2004_parameters = OrderedDict([
            ('rc',2.5),
            ('hc',0.5),
            ('h0',0.5)
        ])
    r = create_r(rmax=r_max,n=r_N)
    dens_w_cutoff = func_cutoff_mishin2004(r=r,**cutoff_mishin2004_parameters)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,1)
    ax.plot(r,dens_w_cutoff)
    plt.show()
