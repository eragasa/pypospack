import unittest
import pypospack.io.vasp as vasp


incar_fin = 'INCAR.in'
incar_fout = 'INCAR.out'

incar = vasp.Incar()
incar.read(incar_fin)

incar.write(incar_fout)

incar.read(incar_fout)
