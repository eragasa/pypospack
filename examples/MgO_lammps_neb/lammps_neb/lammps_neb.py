import sys
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
                import pypospack.vasp import Poscar
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

def file_io__write_header(parameter_filename,filename,qoiname_neb):
    assert type(filename) is str
    assert type(parameter_names) is list
    assert all([type(p) is str for p in parameter_names])
    assert qoiname_neb

    lines=None
    with open(parameter_filename,'r') as f:
        lines = f.readlines()

    names=[str(n.strip()) for n in lines[0].strip().split(',')]
    types=[str(t.strip()) for t in lines[1].strip().split(',')]
    parameter_names=[names[i] for i,v in enumerate(types) if v == 'param']
    qoi_names=[names[i] for i,v in enumerate(types) if v=='qoi']

    names = ['sim_id'] + parameter_names + [qoiname_neb]
    types = ['sim_id'] + len(parameter_names)*['param'] + ['qoi']

    _strout = ",".join(names) + "\n"
    _strout += ",".join(types) + "\n"

    with open(filename,'w') as f:
        f.write(_strout)

def file_io__write_result(filename,sim_id,parameters,neb_value):
    assert type(filename) is str
    assert type(parameters) is OrderedDict
    assert type(neb_value) is float

    values = []
    values.append(sim_id)
    for p in parameters:
        values.append(parameters[p])
    values.append(neb_value)

    _strout = ",".join([str(v) for v in values]) + "\n"

    with open(filename,'a') as f:
        f.write(_strout)

def process_neb_line():
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

    n_image = len(pe)
    e_0 = pe[0]
    e_f = pe[n_image - 1]
    e_max = max(pe)
    e_barrier = e_max - 0.5*(e_0 + e_f)
    return (e_0,e_f,e_max,e_barrier)

if __name__ == "__main__":
    symbols=['Mg','O']
    rcut=10.0
    neb_state=sys.argv[1]
    i_sim=int(sys.argv[2]

    output_dir='output'
    lmps_script_dir='input'
    structure_dir='input'
    lmps_structure_bulk='MgO_NaCl_333.structure'
    lmps_structure_init='MgO_NaCl_333_fr_a_0.structure'
    lmps_structure_final='MgO_NaCl_333_fr_a_1.structure'
    pypospack_param=os.path.join('input','subselect.d_metric.sub_b_lt_median.out')
    pypospack_out=os.path.join('output',
    qoiname_neb = 'MgO_fr_a.neb'
    if neb_state=='start':
        import os,shutil

        # create output directory
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        
        # create lammps structure file: bulk unit cell
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_bulk),
            dst=os.path.join(os.cwd(),lmps_structure_bulk))
        # create lammps structue file: neb initial image
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_init),
            dst=os.path.join(os.cwd(),lmps_structure_init))
        # create lammps structure file: neb final image
        shutil.copyfile(
            src=os.path.join(structure_dir,lmps_structure_final),
            dst=os.path.join(os.cwd(),lmps_structure_final))

        # create lammps script, initial minimization: lmps_min
        shutil.copyfile(
            src=os.path.join(lmps_structure_dir,'in.min'),
            dst=os.path.join(os.cwd(),'in.min'))
        # create lammps_script, neb calculation: in.neb
        shutil.copyfile(
            src=os.path.join(lmps_structure_dir,'in.neb'),
            dst=os.path.join(os.cwd(),'in.neb'))

        # write header to the output file
        file_io__write_header(
            filename_out=pypospack_out,
            filename_param=pypospack_param,
            qoiname_neb=qoi_name) 

    elif neb_state=='next_param':

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
          
    elif neb_stat=='post':
        (e_0,e_f,e_max,e_barrier)=process_neb_line()
        
        file_io__write_result(
                filename_out,
                filename_params
                qoiname_neb=qoiname_neb,
                neb_value=e_barrier)
