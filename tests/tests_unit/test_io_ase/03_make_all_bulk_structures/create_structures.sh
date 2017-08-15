#!/bin/bash
rm *.vasp
python Ni_fcc_unit.py
python Ni_bcc_unit.py
python Ni_hcp_ortho.py
python Ni_dia_unit.py
python Ni_sc_unit.py
