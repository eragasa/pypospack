import sys, os, shutil, re
from collections import OrderedDict

#if neb_state in ['start']:
#    pypospack_parameter_filename=sys.argv[2]
#
#if neb_state in ['next_param']:
#    lmps_structure_bulk=sys.argv[2]
#    lmps_structure_init=sys.argv[3]
#    lmps_structure_final=sys.argv[4]
#    pypospack_parameter_filename=sys.argv[5]
#    pypospack_output_filename='pypospack.neb.out'
#    isim=int(sys.argv[6])

class WorkflowLammpsNebSimulation(object):
    def __init__(self,structure_dir,structures):
        self.structure_dir=structure_dir
        self.structures = OrderedDict()
        self.structures['ideal']=lmps_structure_bulk
        self.structures['neb_initial']=lmps_structure_init
        self.structures['neb_final']=lmps_structure_final

    def write_lammps_structure_files(self):
        for k,v in self.structures:
            if k.endswith('.vasp'):
                #import pypospack.vasp import Poscar
                poscar = Poscar()
                poscar.open(k)

            elif k.endswith('.structure'):
                pass
            else:
                raise ValueError("Unknown structure type")

def get_parameter_set(parameter_filename,i_sim):
    
    lines=None
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

    return parameters

def file_io__write_header(filename_out,filename_param,qoiname_neb):
    assert type(filename_out) is str
    assert type(filename_param) is str
    assert type(qoiname_neb) is str

    # read in parameter file
    lines=None
    with open(filename_param,'r') as f:
        lines = f.readlines()

    # process the header lines: lines[0] and lines[1]
    names=[str(n.strip()) for n in lines[0].strip().split(',')]
    types=[str(t.strip()) for t in lines[1].strip().split(',')]

    # process the parameter names
    parameter_names=[names[i] for i,v in enumerate(types) if v == 'param']
    
    # process the qoi names
    qoi_names=[names[i] for i,v in enumerate(types) if v=='qoi']

    # <------------ create header lines for the output file
    _names_out = ['sim_id'] + parameter_names + [qoiname_neb]
    _types_out = ['sim_id'] + len(parameter_names)*['param'] + ['qoi']

    _strout = ",".join(_names_out) + "\n"
    _strout += ",".join(_types_out) + "\n"

    print('pypospack outfile:{}'.format(filename_out))
    print('header_lines:\n{}'.format(_strout))
    
    with open(filename_out,'w') as f:
        f.write(_strout)

def file_io__write_result(
        filename_out,
        filename_params,
        qoiname_neb,
        neb_value,
        i_sim):
    print('filename_out:',filename_out)
    print('filename_params:',filename_params)
    print('qoiname_neb',qoiname_neb)
    print('neb_value:',neb_value,type(neb_value))
    print('i_sim',i_sim)
    #assert type(filename_out) is str
    #assert type(filename_params) is str
    #assert type(qoiname_neb is str
    #assert type(neb_value) is float
    #assert type(i_sim) is int

    lines = None
    with open(filename_params,'r') as f:
        lines = f.readlines()

    names=[str(n.strip()) for n in lines[0].strip().split(',')]
    types=[str(t.strip()) for t in lines[1].strip().split(',')]
    
    parameter_names=[names[i] for i,v in enumerate(types) if v == 'param']
    qoi_names=[names[i] for i,v in enumerate(types) if v=='qoi']
    values=[n.strip() for n in lines[2+i_sim].strip().split(',')]

    sim_id=values[names.index('sim_id')]
    parameters=OrderedDict(
            [(p,values[names.index(p)]) for p in parameter_names])
    qoi_values=[neb_value]

    if type(sim_id) is str:
        _values = [sim_id]\
                + [str(v) for k,v in parameters.items()]\
                + [str(v) for v in qoi_values]
    else:
        _values = [str(sim_id)] \
                + [str(v) for k,v in parameters.items()]\
                + [str(v) for v in qoi_values]
    
    _strout = ",".join(_values) + "\n"

    with open(filename_out,'a') as f:
        f.write(_strout)

def process_neb_line(lmps_log_file='log.lammps'):
    import re

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

    n_image = len(pe)
    e_0 = pe[0]
    e_f = pe[n_image - 1]
    e_max = max(pe)
    e_barrier = e_max - 0.5*(e_0 + e_f)
    return (e_0,e_f,e_max,e_barrier)

if __name__ == "__main__":
    import os,shutil
    symbols=['Mg','O']
    rcut=10.0

    output_dir='output'
    lmps_script_dir='input'
    structure_dir='input'
    lmps_structure_bulk='MgO_NaCl_333.structure'
    lmps_structure_init='MgO_NaCl_333_fr_a_0.structure'
    lmps_structure_final='MgO_NaCl_333_fr_a_1.structure'
    pypospack_param=os.path.join(
            'input',
            'subselect.d_metric.sum_b_lt_median.out')
    pypospack_out=os.path.join(
            'output',
            'pypospack.neb.out')
    qoiname_neb='MgO_fr_a.neb'
    
    neb_state=sys.argv[1]
    if neb_state=='start':

        # create output directory
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        
        # create lammps structure file: bulk unit cell
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_bulk),
            dst=os.path.join(os.getcwd(),lmps_structure_bulk))
        # create lammps structue file: neb initial image
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_init),
            dst=os.path.join(os.getcwd(),lmps_structure_init))
        # create lammps structure file: neb final image
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_final),
            dst=os.path.join(os.getcwd(),lmps_structure_final))

        # create lammps script, initial minimization: lmps_min
        shutil.copyfile(
            src=os.path.join(structure_dir,'in.min'),
            dst=os.path.join(os.getcwd(),'in.min'))
        # create lammps_script, neb calculation: in.neb
        shutil.copyfile(
            src=os.path.join(structure_dir,'in.neb'),
            dst=os.path.join(os.getcwd(),'in.neb'))

        # write header to the output file
        file_io__write_header(
            filename_out=pypospack_out,
            filename_param=pypospack_param,
            qoiname_neb=qoiname_neb) 

    elif neb_state=='next_param':
        i_sim=int(sys.argv[2])

        # write potential.mod
        parameters=get_parameter_set(pypospack_param,i_sim)
        # potential section
        from pypospack.potential import BuckinghamPotential
        buck = BuckinghamPotential(symbols=symbols)
        _strout = buck.lammps_potential_section_to_string(parameters=parameters,rcut=rcut)
        
        # charge summation method
        _strout += "\n"
        _strout += "kspace_style pppm 1.5e-5\n"
        _strout += "\n"
        _strout += "neighbor 1.0 bin\n"
        _strout += "neigh_modify every 1 delay 0 check yes"
        _fname_out = 'potential.mod'
        with open(_fname_out,'w') as f:
            f.write(_strout)
          
    elif neb_state=='post':
        i_sim=int(sys.argv[2])
        (e_0,e_f,e_max,e_barrier)=process_neb_line()
        print(e_0,e_f,e_max,e_barrier) 
        file_io__write_result(
                filename_out=pypospack_out,
                filename_params=pypospack_param,
                qoiname_neb=qoiname_neb,
                neb_value=e_barrier,
                i_sim=i_sim)
