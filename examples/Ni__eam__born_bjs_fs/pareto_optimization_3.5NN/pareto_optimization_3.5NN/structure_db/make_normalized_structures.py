import argparse
from pypospack.io.vasp import Poscar
_description = 'Normalizes the H-matrix'
_parser = argparse.ArgumentParser(description=_description)
_parser.add_argument('--in', action='store', dest='filename_in', type=str)
_parser.add_argument('--out', action='store', dest='filename_out', type=str)

_args = _parser.parse_args()
_filename_in = _args.filename_in
_filename_out = _args.filename_out

print('filename_in:{}'.format(_filename_in))
_poscar = Poscar()
_poscar.read(_filename_in)
print('a0:{}'.format(_poscar.a0))
print('H:\n{}'.format(_poscar.H))
print('filename_out:{}'.format(_filename_out))
_poscar.normalize_h_matrix()
print(80*'-')
print('a0:{}'.format(_poscar.a0))
print('H:\n{}'.format(_poscar.H))
_poscar.write(_filename_out)

