import os,shutil

from ase import Atom,Atoms
from ase.lattice.cubic import FaceCenteredCubic
import ase.build.bulk

from pypospack.crystal import SimulationCell
from pypospack.io.vasp import Poscar
from pypospack.io.lammps import LammpsStructure
import pypospack.crystal as crystal

def make_Si_sc_bulk(name='Si'):
    name = 'Si'
    cubic = True
    return ase.build.bulk(name,cubic=True)


def make_fcc_bulk(name='Cu'):
    return ase.build.bulk(name,cubic=True)

def make_fcc_111_slab(size=(2,2,3),symbol='Cu',pbc=(1,1,0)):
    direction_x = [1,-1,1]
    direction_y = [1,1,-2]
    direction_z = [1,1,1]
    atoms = FaceCenteredCubic(\
            directions=[direction_x,
                        direction_y,
                        direction_z],
            size=size,
            symbol=symbol,
            pbc=pbc)
    return atoms

def make_fcc_100_slab(size=(2,2,3),symbol='Cu',pbc=(1,1,0)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]

    atoms = FaceCenteredCubic(\
            directions=[direction_x,
                        direction_y,
                        direction_z],
            size=size,
            symbol=symbol,
            pbc=pbc)

if __name__ == "__main__":
    from collections import OrderedDict
    
    is_write_lammps_structures = True
    is_write_vasp_structures = True 
    dst_structure_dir = 'Si_structures'
    structures = OrderedDict()
    structures['Si_dia_unit'] = OrderedDict()
    structures['Si_dia_unit']['ase_obj'] = make_Si_sc_bulk()
    structures['Si_dia_unit']['symbol_list'] = ['Si']
    structures['Si_dia_unit']['lmps_atom_style'] = 'atomic'


    if os.path.exists(dst_structure_dir):
        print('dst_structure_dir exists...')
        print('...deleting directory and contents {}'.format(dst_structure_dir))
        print('...recreating the directory {}'.format(dst_structure_dir))
        shutil.rmtree(dst_structure_dir)
        os.mkdir(dst_structure_dir)
    else:
        print('dst_structure_dir does not exist...')
        print('...creating the directory {}'.format(dst_structure_dir))
        os.mkdir(dst_structure_dir)

    for k,v in structures.items():
        # convert ase -> pypospack
        v['pypospack_obj'] = SimulationCell(v['ase_obj'])
        # convert pypospack -> VASP poscar
        if is_write_vasp_structures:
            v['vasp_obj'] = Poscar(v['pypospack_obj'])
            v['vasp_fname'] = '{}.vasp'.format(k)
            v['vasp_obj'].write(os.path.join(
                dst_structure_dir,
                v['vasp_fname']))
            print('{}->{}'.format(k,v['vasp_fname'])) 
        
        # convert pypospack -> LAMMPS structure
        if is_write_lammps_structures:
            v['lammps_obj'] = LammpsStructure(v['pypospack_obj'])
            v['lammps_fname'] = '{}.structure'.format(k)
            v['lammps_obj'].write(
                filename=os.path.join(
                    dst_structure_dir,
                    v['lammps_fname']
                ),
                symbol_list=v['symbol_list'],
                atom_style=v['lmps_atom_style']
            )
            print('{}->{}'.format(k,v['lammps_fname'])) 

    structures['Si_dia_333'] = OrderedDict()
    structures['Si_dia_333']['pypospack_obj'] = crystal.make_super_cell(
            structure=structures['Si_dia_unit']['pypospack_obj'],
            sc=[3,3,3])
    poscar = Poscar(structures['Si_dia_333']['pypospack_obj'])
    poscar.write(os.path.join(
        dst_structure_dir,
        'Si_dia_333.vasp'))
    lammps_structure = LammpsStructure(structures['Si_dia_333']['pypospack_obj'])
    lammps_structure.write(
            filename=os.path.join(
                dst_structure_dir,
                'Si_dia_333.structure'),
            symbol_list=['Si'],
            atom_style='atomic')

    structures['Si_dia_333_vac'] = OrderedDict()
    structures['Si_dia_333_vac']['pypospack_obj'] = SimulationCell(
            structures['Si_dia_333']['pypospack_obj'])
    structures['Si_dia_333_vac']['pypospack_obj'].add_vacancy(
            symbol = 'Si',
            position = [0.0000,0.0000,0.0000])
    poscar = Poscar(structures['Si_dia_333_vac']['pypospack_obj'])
    poscar.write(os.path.join(
        dst_structure_dir,
        'Si_dia_333_vac.vasp'))
    lammps_structure = LammpsStructure(structures['Si_dia_333']['pypospack_obj'])
    lammps_structure.write(
            filename=os.path.join(
                dst_structure_dir,
                'Si_dia_333_vac.structure'),
            symbol_list=['Si'],
            atom_style='atomic')

