from collections import OrderedDict
import matplotlib.pyplot as plt
import pyflamestk.potential as potential
import pyflamestk.eam as eam
import numpy as np

# configuration inforamtion
symbols = ['Ni']

#<--- variables unique for the test ---------------------------------------
symbols = ['Ni']
func_pair='morse'
func_density='eam_dens_exp'
func_embedding='eam_embed_universal'

#<--- setup of the code to conduct the test -------------------------------
from pypospack.potential import EamPotential

#<--- code being tested ---------------------------------------------------
eam_pot = EamPotential(
        symbols=symbols,
        func_pair=func_pair,
        func_density=func_density,
        func_embedding=func_embedding)

lattice_info = OrderedDict()
lattice_info['Ni'] = OrderedDict()
lattice_info['Ni']['a0'] = 3.524 # in Angs
lattice_info['Ni']['lattice_type'] = 'fcc'

# create parameter diction
parameters = OrderedDict()
parameters['p_NiNi_D0'] = 2.
parameters['p_NiNi_a'] = 1
parameters['p_NiNi_r0'] = 1
parameters['d_Ni_rho0'] = 1
parameters['d_Ni_beta'] = 1
parameters['d_Ni_r0'] = 1
parameters['e_Ni_F0'] = 1
parameters['e_Ni_p'] = 1/2   # this is equilvalent to FS 
parameters['e_Ni_q'] = 1.00   # this is equilvalent to FS
parameters['e_Ni_F1'] = 1
parameters['e_Ni_rho0'] = 1
r_cutoffs = OrderedDict()

fcc_NN_info = OrderedDict()
fcc_NN_info['1NN'] = OrderedDict()
fcc_NN_info['2NN'] = OrderedDict()
fcc_NN_info['3NN'] = OrderedDict()
fcc_NN_info['1NN']['n'] = 12
fcc_NN_info['2NN']['n'] = 6
fcc_NN_info['3NN']['n'] = 24
fcc_NN_info['1NN']['da'] = (2**0.5)/2
fcc_NN_info['2NN']['da'] = 1.
fcc_NN_info['3NN']['da'] = 1.225

Nr = 5000
Nrho = 5000

for s in symbols:
    if lattice_info[s]['lattice_type'] == 'fcc':

        a0 = lattice_info[s]['a0']
        da_2NN = fcc_NN_info['2NN']['da']
        da_3NN = fcc_NN_info['3NN']['da']
        r_cutoffs[s] = 0.5 * a0 * (da_2NN+da_3NN)
print(80*'-')
print('{:^80}'.format('RADIAL CUTOFF'))
print(80*'-')
for s in symbols:
    print('{:>2} {:^8} {:+10.6e}'.format(s,lattice_info[s]['lattice_type'],r_cutoffs[s]))    

r0 = OrderedDict()
for s in symbols:
    if lattice_info[s]['lattice_type'] == 'fcc':
        r0[s] = fcc_NN_info['1NN']['da'] * lattice_info[s]['a0']
parameters['p_NiNi_r0']= r0['Ni']
parameters['d_Ni_r0'] = r0['Ni']

#print(80*'-')
#print('{:^80}'.format('RADIUS VECTOR'))
#print(80*'-')
rmax = max([v for k,v in r_cutoffs.items()])
r = rmax * np.linspace(1,Nr,Nr)/Nr
dr = r[1]-r[0]


eam_pot.evaluate_density(r=r,parameters=parameters,rcut=r_cutoffs['Ni'])
#plt.plot(r,eam_pot.density['Ni'])
#plt.savefig('density.png')

# determine equilibrium density
print(80*'-')
print('{:^80}'.format('DETERMINE EQUILIBRIUM DENSITY'))
print(80*'-')
rho0 = OrderedDict()
for s in symbols:
    print('elements:{}'.format(s))
    rho0[s] = 0
    for iNN in ['1NN','2NN']:
        n = fcc_NN_info[iNN]['n']
        d = fcc_NN_info[iNN]['da']*lattice_info['Ni']['a0']
        d_rho0 = np.interp(x=d,xp=r,fp=eam_pot.density['Ni'])
        rho0[s] += n*d_rho0
        print('{} {:>2} {:+10.6e} {:+10.6e} {:+10.6e}'.format(
            iNN, n,d,d_rho0,rho0[s]))
    print(80*'-')
for s in symbols:
    print('rho0({})={:+10.6e}'.format(s,rho0[s]))

print(80*'-')
print('{:^80}'.format('DETERMINE MAXIMUM DENSITY'))
print(80*'-')
max_compression = 0.50
print('max_compression = {}'.format(max_compression))
rhomax = OrderedDict()
for s in symbols:
    print('elements:{}'.format(s))
    rhomax[s] = 0
    for iNN in ['1NN','2NN']:
        n = fcc_NN_info[iNN]['n']
        d = (1-max_compression)*fcc_NN_info[iNN]['da']*lattice_info['Ni']['a0']
        d_rho0 = np.interp(x=d,xp=r,fp=eam_pot.density['Ni'])
        rhomax[s] += n*d_rho0
        print('{} {:>2} {:+10.6e} {:+10.6e} {:+10.6e}'.format(
            iNN, n,d,d_rho0,rho0[s]))
    print(80*'-')
for s in symbols:
    print('rhomax({})={:+10.6e}'.format(s,rhomax[s]))
rhomax = max([v for k,v in rhomax.items()])
print('rhomax={:10.6e}'.format(rhomax))

# change the parameters
s = 'Ni'
print(80*'-')
print('{:^80}'.format('PARAMETERS'))
print(80*'-')
parameters['e_{}_rho0'.format(s)] = rho0[s]
for k,v in parameters.items():
    print('{:^20} {:+10.6e}'.format(k,v))
# make the equilibrium value twice the value of the fcc lattice

rho = rhomax * np.linspace(1,Nrho,Nrho)/Nrho
eam_pot.evaluate_embedding(rho=rho,parameters=parameters)
embedding_fs = 2*parameters['e_Ni_F0']*np.sqrt(rho/parameters['e_Ni_rho0']) 

#plt.plot(rho,eam_pot.embedding['Ni'])
#plt.plot(rho,embedding_fs)
#plt.savefig('eam_embedding.png')
#exit()

fcc_a0= lattice_info['Ni']['a0']
fcc_a_min = 0.80 * fcc_a0
fcc_a_max = 1.25 * fcc_a0
#fcc_a_min = 5
#fcc_a_max = 6
fcc_a = fcc_a_min + np.linspace(1,Nr,Nr)/Nr*(fcc_a_max-fcc_a_min)
print('fcc_a_min:{},fcc_a_max:{}'.format(fcc_a.min(),fcc_a.max()))
print('fcc_a: {} {}'.format(fcc_a.min(),fcc_a.max()))

# fcc_rho:  rho(fcc_r)
fcc_rho = []
for _a in fcc_a:
    _rho = 0
    _phi = 0
    for iNN in ['1NN','2NN']: 
        n = fcc_NN_info[iNN]['n']
        d = fcc_NN_info[iNN]['da']*_a
        d_rho0 = np.interp(x=d,xp=r,fp=eam_pot.density['Ni'])
        _rho += n*d_rho0
    fcc_rho.append(_rho)
fcc_rho = np.array(fcc_rho)

#plt.plot(fcc_a,fcc_rho)
#plt.savefig('dens_fcc_a.png')
#exit()
import copy
fcc_embed = np.interp(x=fcc_rho,xp=rho,fp=embedding_fs)
eam_pot.evaluate_pair(r=r,parameters=parameters,rcut=r_cutoffs['Ni'])
fcc_phi = []
for _a in fcc_a:
    _phi = 0
    for iNN in ['1NN','2NN']:
        n = fcc_NN_info[iNN]['n']
        d = fcc_NN_info[iNN]['da']*_a
        d_phi0 = np.interp(x=d,xp=r,fp=eam_pot.pair['NiNi'])
        _phi += n*d_phi0
    fcc_phi.append(_phi)
fcc_phi = np.array(fcc_phi)

f,ax = plt.subplots(2)
phi_cut = copy.deepcopy(eam_pot.pair['NiNi'])
eam_pot.evaluate_pair(r,parameters)
phi_nocut = copy.deepcopy(eam_pot.pair['NiNi'])
ax[0].plot(r,phi_cut)
ax[0].plot(r,phi_nocut)
ax[0].set_xlim(2.,r.max())
ax[0].set_ylim( 
        min(phi_cut.min(),phi_nocut.min()),
        min(np.interp(2.,r,phi_cut),np.interp(2.,r,phi_nocut))
    )
ax[1].plot(fcc_a,fcc_phi)
ax[1].plot(fcc_a,fcc_embed)
ax[1].plot(fcc_a,fcc_phi+fcc_embed)
f.savefig('eam_fcc.png')
exit()
#plt.plot(fcc_a,fcc_phi)
#plt.savefig('fcc_pair.png')
#exit()

#plt.plot(fcc_a,fcc_embed)
#plt.xlabel('fcc a')
#plt.ylabel('embedding energy')
#plt.savefig('embed_fcc_a.png')
#exit()
#fcc_embed = np.interp(x=fcc_rho,xp=rho,fp=eam_pot.embedding['Ni'])

print('fcc_a_type:{}'.format(str(type(fcc_a))))
print('fcc_rho_type:{}'.format(str(type(fcc_rho))))
print('rho: min={}, max={}'.format(rho.min(),rho.max()))
print('fcc_rho: min={}, max={}'.format(fcc_rho.min(),fcc_rho.max()))
print('fcc_embed: min={}, max={}'.format(fcc_embed.min(),fcc_embed.max()))
#plt.plot(fcc_a,fcc_embed)
#plt.ylim(fcc_embed.min(),fcc_embed.max())
#plt.savefig('eam_embedding.png')
#exit()
#plt.plot(r,eam_pot.density['Ni'])
#d_parameter_names = []
#for s in symbols:
#    d_parameter_names += eam_pot.obj_density.parameter_names
#d_parameters = OrderedDict()
#for d_pn in d_parameter_names:
#    pn = 'd_{}'.format(d_pn)
#    d_parameters[pn] = parameters[pn]

# generate evently spaced NRho points between rho_low, and rho_high
# units for rho is unclear, but should be a measure of electron density
rho_low = 0
rho_high = rhocut_g
rho = np.linspace(rho_low,rho_high,Nrho)
drho = rho[1]-rho[0]

dens_parameters = {}
embed_parameters = {}
for s in symbols:
    dens_parameters[s] = {}
    embed_parameters[s] = {}
    for k,v in parameters.items():
        if k.startswith('e.{}'.format(s)):
            key = k.replace('e.','')
            embed_parameters[s][key] = v
        elif k.startswith('d.{}'.format(s)):
            key = k.replace('d.','')
            dens_parameters[s][key] = v

pair_parameters = {}
N_symbols = len(symbols)
for i in range(N_symbols):
    for j in range(N_symbols):
        if i >= j:
            pair = '{}{}'.format(symbols[i],symbols[j])
            pair_parameters[pair] = {}
            for k,v in parameters.items():
                if k.startswith('p.{}'.format(pair)):
                    key = k.replace('p.','')
                    pair_parameters[pair][key] = v

print(dens_parameters)
print(embed_parameters)
print(pair_parameters)

