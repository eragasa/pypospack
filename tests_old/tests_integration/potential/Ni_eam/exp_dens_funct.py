import numpy as np
import matplotlib.pyplot as plt
import pyflamestk.potential as potential

# ---- test for one symbol
symbols = ['Ni']
expfunc = potential.ExponentialDensityFunction(symbols)
print('--- test exponential density function for one symbol')
print(expfunc.potential_type)
print(expfunc.parameter_names)

# --- test for two symbols
symbols = ['Ni','Al']
expfunc = potential.ExponentialDensityFunction(symbols)
print('--- test exponential density function for two symbols')
print(expfunc.potential_type)
print(expfunc.parameter_names)

# --- parameters for test
symbols = ['Ni']
param_dict = {}
param_dict['Ni_rho0'] = 1
param_dict['Ni_beta'] = 1
param_dict['Ni_r0'] = 1

r_low = 0.5
r_high = 10
N_r = 1000
rcut = 4
h = 1
r = np.linspace(r_low, r_high, N_r)

# ---- test with no cutoff function
print('--- calculate without cutoff')
symbol = 'Ni'
expfunc_nocut = potential.ExponentialDensityFunction(symbols)
vals_nocut = expfunc_nocut.evaluate(r=r,symbol='Ni',params=param_dict)

print('--- calculate with cutoff')
symbol = 'Ni'
expfunc_cut = potential.ExponentialDensityFunction(symbols)
vals_cut = expfunc_cut.evaluate(r=r,symbol='Ni',params=param_dict,
        rcut=rcut,h=h)

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(r,vals_nocut,color='b',label='no_cut')
ax.plot(r,vals_cut,color='r',label='cut')
fig.savefig('exp_dens_funct.png',dpi=800)
plt.close(fig)
