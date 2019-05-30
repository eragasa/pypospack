import numpy as np
rmax=15
npoints=1000
r = rmax * np.linspace(1,npoints,npoints)/npoints
dr_step = r[1] - r[0]
dr_theo = rmax/(npoints)
print(r)
print('r.shape={}'.format(r.shape))
print('dr_step={}'.format(dr_step))
print('dr_theo={}'.format(dr_theo))

self.V = self.V - (self.V[rcut_idx] - self.V[cut_idx-1])/self.dr)\
                  *(self.r-self.r_cut)
rcut=10
rcut_idx = max(np.where(r <= rcut)[0])
print(rcut_idx)
print(type(rcut_idx))
