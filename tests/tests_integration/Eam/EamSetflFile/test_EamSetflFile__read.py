import pytest
import os
import numpy as np
from collections import OrderedDict
import pypospack.eamtools as eamtools

def almost_equal(x1,x2,threshold=6):
    episilon = sum(x1,x2)/2 * float("1e-{}".format(threshold))
    return abs(x1-x2) < epsilon

class Test__read__not_auto(object):

    def setup__setup_variables(self):
        self.filename = '../Ni1_Mendelev_2010.eam.fs.txt'
        self.symbols = ['Ni']
        self.n_symbols = 1
        self.N_rho = 10000
        self.N_r = 10000
        self.d_rho = 0.02
        self.d_r = 0.0006
    
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

        assert self.eamfile.N_VALUES_PER_LINE_RHO == 5
        assert self.eamfile.N_VALUES_PER_LINE_R == 5
        assert self.eamfile.SETFL_NUM_FORMAT == "{:+24.16E}"
        assert self.eamfile.SETFL_INT_FORMAT == "{:5d}"
        assert self.eamfile.PAIR_KEY_FORMAT == "{}.{}"
        assert self.eamfile.N_headersection_lines == 5
        assert self.eamfile.N_density_lines is None
        assert self.eamfile.N_embedding_lines is None
        assert self.eamfile.N_pairpotential_lines is None
        assert self.eamfile.N_atomicsection_lines is None
        assert self.eamfile.N_pairpotential_lines is None
        assert isinstance(self.eamfile.comments,list)
        assert self.eamfile.filename is None
        assert self.eamfile.symbols is None
        assert self.eamfile.n_symbols is None
        assert self.eamfile.symbol_dict is None
        assert self.eamfile.r is None
        assert self.eamfile.rho is None
        assert self.eamfile.func_embedding is None
        assert self.eamfile.func_pairpotential is None
        assert self.eamfile.func_density is None
        assert self.eamfile.N_rho is None
        assert self.eamfile.d_rho is None
        assert self.eamfile.max_rho is None
        assert self.eamfile.N_r is None
        assert self.eamfile.d_r is None
        assert self.eamfile.max_r is None
        assert self.eamfile.rcut_g is None

    def test__read_file(self):
        self.setup__initialize()
        self.eamfile = eamtools.EamSetflFile()
        self.eamfile.read(
                filename=self.filename,
                autoprocess=False)

        assert os.path.abspath(self.eamfile.filename) \
                == os.path.abspath(self.filename)
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

        assert self.eamfile.symbols == self.symbols
        assert self.eamfile.n_symbols == self.n_symbols

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
    
        assert type(self.eamfile.N_rho) is int
        assert type(self.eamfile.N_r) is int
        assert type(self.eamfile.d_rho) is float
        assert type(self.eamfile.d_r) is float
        assert type(self.eamfile.r_cut) is float

        assert self.eamfile.symbols == ['Ni']
        assert self.eamfile.n_symbols == 1
        assert self.eamfile.N_rho == 10000
        assert self.eamfile.N_r == 10000
        assert self.eamfile.d_rho == 0.02
        assert self.eamfile.d_r == 0.0006

    def test__process_setfl_header_section(self):
        self.setup__setup_variables()
        eamfile = eamtools.EamSetflFile()
        eamfile.read(filename=self.filename,autoprocess=False)
        eamfile.process_setfl_header_section()
        
        assert type(eamfile.N_rho) is int
        assert type(eamfile.N_r) is int
        assert type(eamfile.d_rho) is float
        assert type(eamfile.d_r) is float
        assert type(eamfile.r_cut) is float

        assert eamfile.symbols == ['Ni']
        assert eamfile.n_symbols == 1
        assert eamfile.N_rho == 10000
        assert eamfile.N_r == 10000
        assert eamfile.d_rho == 0.02
        assert eamfile.d_r == 0.0006
    
        assert eamfile.N_headersection_lines == 5
        assert eamfile.N_density_lines == 2000
        assert eamfile.N_embedding_lines == 2000
        assert eamfile.N_pairpotential_lines == 2000
        assert eamfile.N_atomicsection_lines == 4001
        assert eamfile.N_pairsection_lines == 2000
        n_total_lines = eamfile.N_headersection_lines\
                + eamfile.N_atomicsection_lines\
                + eamfile.N_pairsection_lines
        assert n_total_lines == 6006
        assert len(eamfile.lines) == n_total_lines

    def test__process_setfl_atomic_section(self):
        self.setup__setup_variables()
        eamfile = eamtools.EamSetflFile()
        eamfile.read(filename=self.filename,autoprocess=False)
        eamfile.process_setfl_header_section()
        eamfile.process_setfl_atomic_section()
        
        assert isinstance(eamfile.func_embedding,dict)
        for s in eamfile.symbols:
            assert isinstance(eamfile.func_embedding[s],np.ndarray)
            (n_row,n_col) = np.atleast_2d(eamfile.func_embedding[s]).T.shape
            assert eamfile.N_rho == n_row
        
        assert isinstance(eamfile.func_density,dict)
        for s in eamfile.symbols:
            assert isinstance(eamfile.func_embedding[s],np.ndarray)
            (n_row,n_col) = np.atleast_2d(eamfile.func_density[s]).T.shape
            assert eamfile.N_r == n_row

        assert eamfile.symbol_dict['Ni']['atomic_number'] == 28
        assert eamfile.symbol_dict['Ni']['atomic_mass'] == 58.71
        assert eamfile.symbol_dict['Ni']['lattice_constant'] == 3.518121
        assert eamfile.symbol_dict['Ni']['lattice_type'] == 'fcc'

    def test__process_setfl_pairpotential_section(self): 
        self.setup__setup_variables()
        eamfile = eamtools.EamSetflFile()
        eamfile.read(filename=self.filename,autoprocess=False)
        eamfile.process_setfl_header_section()
        eamfile.process_setfl_atomic_section()
        eamfile.process_setfl_pairpotential_section()

        assert isinstance(eamfile.func_pairpotential,OrderedDict)
        for pair,pairpot in eamfile.func_pairpotential.items():
            assert isinstance(pairpot,np.ndarray)
            (n_row,n_col) = np.atleast_2d(pairpot).T.shape
            assert eamfile.N_r == n_row

class Test__read__auto(object):

    def setup__setup_variables(self):
        self.filename = '../Ni1_Mendelev_2010.eam.fs.txt'
        self.symbols = ['Ni']
        self.n_symbols = 1
        self.N_rho = 10000
        self.N_r = 10000
        self.d_rho = 0.02
        self.d_r = 0.0006

    def test__read(self):
        self.setup__setup_variables()
        eamfile = eamtools.EamSetflFile()
        eamfile.read(filename=self.filename,autoprocess=True)
