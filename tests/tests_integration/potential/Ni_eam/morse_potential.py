import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyflamestk.potential as potential

# ----- test for one symbol
symbols = ['Ni']
morse = potential.MorsePotential(symbols)

print(morse.potential_type)
print(morse.parameter_names)

# ------ test for two symbols
symbols = ['Ni','Al']
morse = potential.MorsePotential(symbols)
print(morse.potential_type)
print(morse.parameter_names)

# ------ parameters for test
symbols = ['Ni']
pair = ['Ni','Ni']
param_dict = {}
param_dict['NiNi_D0'] = 2
param_dict['NiNi_a']  = 0.25
param_dict['NiNi_r0'] = 3

# ------ define distances for numpy array
r_low = 0.5 # Angstom
r_high = 10 # Angstrom
N_r = 1000
rcut = 8
h = 1
r = np.linspace(r_low, r_high, N_r)

# ----- test with no cutoff function
print("calculate without cutoff")
morse_nocut = potential.MorsePotential(symbols)
morse_nocut_vals = morse_nocut.evaluate(r=r,pair=pair,params=param_dict)

# ----- test without cutoff function
print("calculate with cutoff")
morse_cut   = potential.MorsePotential(symbols)
morse_cut_vals = morse_cut.evaluate(r=r,pair=pair,params=param_dict,
                                    rcut=rcut,h=h)                      

for n in range(N_r):
    print(n, r[n], morse_nocut_vals[n], morse_cut_vals[n])

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(r,morse_nocut_vals,color = 'b',label='nocut')
ax.plot(r,morse_cut_vals, color='r',label='cut')
ax.set_xlabel(r"interatomic separation $\AA$")
ax.set_ylabel(r"interatomic energy eV")
fig.savefig('morse_potential.png',dpi=800)   # save the figure to file
plt.close(fig)
