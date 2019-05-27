import os
from collections import OrderedDict
from pypospack.eamtools import EamSetflFile


#setfl_filename = os.path.join(
#        'test_EamSetflFile',
#        'Ni.Zhou2004.set')

setfl_filename = os.path.join(
        'test_EamSetflFile',
        'Ni1_Mendelev_2010.eam.fs')
setflfile = EamSetflFile()
setflfile.read(filename=setfl_filename)

print('symbols:{}'.format(setflfile.symbols))
print('N_r:{}'.format(setflfile.N_r))
print('d_r:{}'.format(setflfile.d_r))
print('max_r:{}'.format(setflfile.max_r),setflfile.N_r*setflfile.d_r)
print('N_rho:{}'.format(setflfile.N_rho))
print('d_rho:{}'.format(setflfile.d_rho))
print('max_rho:{}'.format(setflfile.max_rho),setflfile.N_rho*setflfile.d_rho)
print('r_cut:{}'.format(setflfile.r_cut))

r = setflfile.r
rho = setflfile.rho
pair = setflfile.func_pairpotential['Ni.Ni']
embed = setflfile.func_embedding['Ni']
dens = setflfile.func_density['Ni']

print('len(r):{}'.format(len(r)))
print('len(rho):{}'.format(len(rho)))
print('len(pair):{}'.format(len(pair)))
print('len(embed):{}'.format(len(embed)))
print('len(dens):{}'.format(len(dens)))

import matplotlib.pyplot as plt
fig, axarr = plt.subplots(3,1)

r_low = 1
r_high = 10
phi_low = pair.min()
phi_high = 0.5
axarr[0].plot(r,pair/r)
axarr[0].set_xlim(r_low,r_high)
axarr[0].set_ylim(phi_low,phi_high)
axarr[1].plot(rho,embed)
axarr[2].plot(r,dens)
axarr[2].set_xlim(r_low,r_high)
axarr[2].set_ylim(0,10)
fig.savefig('eam.png')
plt.close(fig)
