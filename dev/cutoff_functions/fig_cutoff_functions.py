import numpy as np
import matplotlib.pyplot as plt
from pypospack.potential.pair_lj import func_lj

def func_cutoff(x,y,x_cut):

    # define the cutoff indicator, 1 except when x > x_cut
    cutoff_ind = np.ones(x.size)
    cutoff_ind[x > x_cut] = 0

    assert cutoff_ind.size == y.size
    return  cutoff_ind * y

def create_r_ij(r_max,n):
    r = r_max/n * np.linspace(1,n,n)
    return r

def differentiate(x,y):
    dy = np.diff(y)
    dx = np.diff(x)

    return dy/dx

class TestCutoffFunction(object):

    def __init__(self, r_max=None, r_cut=None, r_N=None):
        self.r_max = r_max
        self.r_cut = r_cut
        self.r_N = r_N

        if self.r_max is not None and self.r__cut is not None:
            self.create_r()

        self.fig = None
        self.ax = None
        self.x_limits =(0.5,3)
        self.y_limits = (-2,2)

    def create_r(self,r_max=None,r_cut=None,r_N=None):

        if r_max is not None: self.r_max = r_max
        if r_cut is not None: self.r_cut = r_cut
        if r_N is not None: self.r_N = r_N

       
        self.r = self.r_max/self.r_N * np.linspace(1,self.r_N,self.r_N)
        return self.r

    def subplots(self):
        plt.close('all')
        self.fig, self.ax = plt.subplots(2,1)

    def plot_potential(self,func,func_parameters):

        x = self.r
        y = func(x,**func_parameters)

        self.ax[0].plot(x,y,alpha=0.5)
        self.ax[0].set_xlim(self.x_limits)
        self.ax[0].set_ylim(self.y_limits)
    
    def plot_forces(self,func,func_parameters):

        x = self.r
        f = func(x,**func_parameters)
        y = differentiate(x,f)
        self.ax[1].plot(x[1:],y, alpha=0.5)

    def plot_potential_w_cutoffs(self,func,func_parameters,r_cut=None):

        if r_cut is not None: self.r_cut = r_cut

        x = self.r
        f = func(x,**func_parameters)
        y = func_cutoff(x,f,r_cut)

        self.ax[0].plot(x,y,alpha=0.5)
        self.ax[0].set_xlim(self.x_limits)
        self.ax[0].set_ylim(self.y_limits)

    def plot_forces_w_cutoffs(self,func,func_parameters,r_cut=None):
        x = self.r
        f = func(x,**func_parameters)
        f = func_cutoff(x,f,r_cut)
        y = differentiate(x,f)

        self.ax[1].plot(x[1:],y,alpha=0.5)
        self.ax[1].set_xlim(self.x_limits)
        self.ax[1].set_ylim(self.y_limits)

    def create_plot(self):

        plt.show()

if __name__ == "__main__":
    r_max = 5.
    r_cut = 2.
    r_N = 1000
    fig_fn = "fig_cutoff_functions.eps"

    o = TestCutoffFunction()
    o.subplots()
    o.create_r(r_max=r_max,r_cut=r_cut,r_N=r_N)
    assert isinstance(r_max,int) or isinstance(r_max,float)
    assert isinstance(r_cut,int) or isinstance(r_max,float)
    assert o.r.size == r_N
    
    lj_parameters = {'epsilon':1,'sigma':1,'r_cut_pair':None}
    o.plot_potential(func=func_lj,func_parameters=lj_parameters)
    o.plot_potential_w_cutoffs(func=func_lj,func_parameters=lj_parameters,r_cut=r_cut)
    o.plot_forces(func=func_lj,func_parameters=lj_parameters)
    o.plot_forces_w_cutoffs(func=func_lj,func_parameters=lj_parameters,r_cut=r_cut)
    plt.show()
    o.fig.savefig(fig_fn,dpi=1200)
    exit()
    plt.close('all')
    r = create_r_ij(r_max,n)
    fig,ax = plt.subplots(2,1)
    
    x_lims=(0.25,r.max())
    y_lims=(-2.,2.)
    ax[0].set_xlim(x_lims)
    ax[0].set_ylim(y_lims)
    ax[1].set_xlim(x_lims)
    ax[1].set_ylim(y_lims)
    
    ax[0].plot(
            x=r,
            y=func_lj(r,**lj_param))
    ax[0].plot(
            x=r,
            y=func_cutoff(r,func_lj(r,**lj_param),r_cut))
    plt.show()
    exit()
    from scipy.misc import derivative
    ax[1].plot(
            x=r,
            y=derivative(
                lambda z: func_lj(z,**lj_param),
                r)
            )
    ax[1].plot(
            x=r,
            y=derivative(
                lambda z: func_cutoff(z,func_lj(z,**lj_param),r_cut),
                r)
            )
    plt.show()
    
    exit()

