import numpy as np
import matplotlib.pyplot as plt
from pypospack.eamtools import create_r
from pypospack.potential.pair_general_lj import func_cutoff_mishin2004

r = create_r(6.,5000)
rc = 5.168
hc = 0.332

xrc = (r-rc)/hc
psirc = (xrc**4)/(1+xrc**4)
rc_ind = np.ones(r.size)
rc_ind[r > rc] = 0

psirc = psirc * rc_ind

h0 = 0.332
x0 = r/h0
psi0 = (x0**4)/(1+x0**4)
fig, ax = plt.subplots(3, 1)

ax[0].plot(r,psirc,label=r'$\Psi_{c}$')
ax[0].set_ylabel(r'$\Psi_{c}$')

ax[1].plot(r,psi0,label=r'$\Psi_{0}$')
ax[1].set_ylabel(r'$\Psi_{0}$')

ax[2].plot(r,psirc*psi0,label=r'$\Psi_{c}*\Psi_{0}$')
ax[2].set_ylabel(r'$\Psi_{c}\Psi_{0}$')

for i in range(2):
    ax[i].tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
fig.tight_layout()
fig.savefig('fig_cutoff_mishin2004.png',dpi=1300)

ax[2].plot(r,
           func_cutoff_mishin2004(r,rc,hc,h0))
plt.show()
