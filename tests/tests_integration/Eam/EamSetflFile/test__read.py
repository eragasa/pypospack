import pytest
import os
import numpy as np

import pypospack.eamtools as eamtools
class Test__read__not_auto(object):

    def setup__setup_variables(self):
        self.filename = '../Ni1_Mendelev_2010.eam.fs.txt'

    def test__testing_harness__setup_variables(self):
        self.setup__setup_variables()
        assert os.path.isfile(self.filename)

    def setup__import_module(self):
        self.setup__setup_variables()
    
    def test__can_import_module(self):
        self.setup__import_module()

    def setup__initialize(self):
        self.setup__import_module()
        self.eamfile = eamtools.EamSetflFile()

    def test__init(self):
        self.setup__initialize()

    def test__read_file(self):
        self.setup__initialize()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)

        assert self.eamfile.filename == self.filename
        assert isinstance(self.eamfile.lines,list)

    def test__process_setfl_comments(self):
        self.setup__setup_variables()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)
        self.eamfile.process_setfl_comments()

        assert type(self.eamfile.comments) is list
        assert len(self.eamfile.comments) ==  3

    def test__process_setfl_n_symbols(self):
        self.setup__setup_variables()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)
        self.eamfile.process_setfl_comments()
        self.eamfile.process_setfl_n_symbols()

        assert type(self.eamfile.symbols) is list
        assert type(self.eamfile.n_symbols) is int
        assert len(self.eamfile.symbols) == self.eamfile.n_symbols
    def test__process_setfl_rho_r_args(self):
        self.setup__setup_variables()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)
        self.eamfile.process_setfl_comments()
        self.eamfile.process_setfl_n_symbols()
        self.eamfile.process_setfl_rho_r_args()

    def test__process_setfl_header_section(self):
        self.setup__setup_variables()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)
        self.eamfile.process_setfl_header_section()

        assert type(self.eamfile.lines) is list
        assert type(self.eamfile.max_rho) is float
        assert type(self.eamfile.max_r) is float
        assert type(self.eamfile.r) is np.ndarray
        assert type(self.eamfile.rho) is np.ndarray

    def test__process_setfl_atomic_section(self):
        self.setup__setup_variables()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)
        self.eamfile.process_setfl_header_section()
        
