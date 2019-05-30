#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import OrderedDict

import numpy as np

import pypospack.utils
from pypospack.eamtools import SeatonSetflReader
from pypospack.potential import EamPotential

setfl_filename = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')

setfl_reader = SeatonSetflReader(path=setfl_filename)
setfl_reader.read()


density_types = ['eam_dens_exp','eam_dens_mishin']
embedding_types = ['eam_embed_universal','eam_embed_bjs'

o_potential = EamPotential(symbols=['Ni','Al'],
        func_pair=)

