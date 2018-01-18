import re

lmps_log_file='log.lammps'

lines=None
with open(lmps_log_file,'r') as f:
    lines=f.readlines()

line = lines[len(lines)-1].strip()
line = re.sub(' +',' ',line)
line = [float(s) for s in line.split(" ")]
n_data = len(line)
step = line[0]
max_replica_force = line[1]
max_atom_force = line[2]
grad_v0 = line[3]
grad_v1 = line[4]
grac_vc = line[5]
ebf = line[6]
ebr = line[7]
rdt = line[8]

rd = []
pe = []
for i in range(9,n_data):
    if i % 2 == 0:
        pe.append(line[i])
    else:
        rd.append(line[i])

#import matplotlib.pyplot as plt
#plt.plot(rd,pe)
#print(rd)
#print(pe)

n_image = len(pe)
e_0 = pe[0]
e_f = pe[n_image - 1]
e_max = max(pe)
e_barrier = e_max - 0.5*(e_0 + e_f)
print(e_0,e_f,e_max,e_barrier)
