import sys,re
from collections import OrderedDict

state=sys.argv[1]
i_sim = int(sys.argv[2])
parameter_filename=sys.argv[3]

class LammpsNebSimulation:
    pass

def get_parameter_set(parameter_filename,i_sym):
    lines = None
    with open(parameter_filename,'r') as f:
        lines = f.readlines()

    names=[str(n.strip()) for n in lines[0].strip().split(',')]
    types=[str(t.strip()) for t in lines[1].strip().split(',')]

    values=[n.strip() for n in lines[2+i_sim].strip().split(',')]
    sim_id = int(float(values[0]))
    param_names = [n for i,n in enumerate(names) if types[i]=='param']
    parameters = OrderedDict()
    for p in param_names:
        v = values[names.index(p)]
        parameters[p] = v

    print(sim_id,parameters)
    return sim_id,parameters

def make_neb_file():
    fname = str(sys.argv[1])
    lines = None
    with open(fname,'r') as f:
        lines = f.readlines()

def process_neb_simulation():

    f = open('log.lammps','r')
    lines = f.readlines()
    f.close()
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

if __name__ == "__main__":
    if state=='min':
        parameters = get_parameter_set(parameter_filename,i_sim)
        print('min')
    elif state=='neb':
        print('neb')
    elif state=='post':
        print('post')
        process_neb_simulation()
    else:
        pass
