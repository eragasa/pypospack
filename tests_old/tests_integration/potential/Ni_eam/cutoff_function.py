import numpy as np
import matplotlib.pyplot as plt
import pyflamestk.potential as potential

# define distances as a numpy array
r_low = 0 # Angstom
r_high = 20 # Angstrom
N_r = 1000
rcut = 10
h = 1
r = np.linspace(r_low, r_high, N_r)

# test cutoff function
phi = {}
for h in [0.01,0.1,1,10]:
    phi[h] = potential.func_cutoff(r,rcut,h)

fig = plt.figure()
ax = plt.subplot(111)
for h in [0.01,0.1,1,10]:
    ax.plot(r,phi[h],label='{}'.format(h))
fig.savefig('cutoff_function.png',dpi=800)
plt.close(fig)


