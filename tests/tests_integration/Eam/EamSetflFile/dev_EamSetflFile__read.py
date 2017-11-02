import pypospack.eamtools as eamtools
import pypospack.io.filesystem as filesystem
from collections import OrderedDict
import numpy as np

if __name__ == "__main__":
    import os

    filename = '../Ni1_Mendelev_2010.eam.fs.txt'
    symbols = ['Ni']
    n_symbols = 1
    assert os.path.isfile(filename)

    # initialization
    eamfile = eamtools.EamSetflFile()
    assert eamfile.N_VALUES_PER_LINE_RHO == 5
    assert eamfile.N_VALUES_PER_LINE_R == 5
    assert eamfile.SETFL_NUM_FORMAT == "{:+24.16E}"
    assert eamfile.SETFL_INT_FORMAT == "{:5d}"
    assert eamfile.PAIR_KEY_FORMAT == "{}.{}"
    assert eamfile.N_headersection_lines == 5
    assert eamfile.N_density_lines is None
    assert eamfile.N_embedding_lines is None
    assert eamfile.N_pairpotential_lines is None
    assert eamfile.N_atomicsection_lines is None
    assert eamfile.N_pairpotential_lines is None
    assert isinstance(eamfile.comments,list)
    assert eamfile.filename is None
    assert eamfile.symbols is None
    assert eamfile.n_symbols is None
    assert eamfile.symbol_dict is None
    assert eamfile.r is None
    assert eamfile.rho is None
    assert eamfile.func_embedding is None
    assert eamfile.func_pairpotential is None
    assert eamfile.func_density is None
    assert eamfile.N_rho is None
    assert eamfile.d_rho is None
    assert eamfile.max_rho is None
    assert eamfile.N_r is None
    assert eamfile.d_r is None
    assert eamfile.max_r is None
    assert eamfile.rcut_g is None

    # read step_though, not autoprocess
    eamfile.read(filename=filename,autoprocess=False)
    assert os.path.abspath(filename) == os.path.abspath(eamfile.filename)
    assert isinstance(eamfile.lines,list)

    # test__process_setfl_comments
    eamfile.process_setfl_comments()
    assert type(eamfile.comments) is list
    assert len(eamfile.comments) == 3
    for c in eamfile.comments:
        assert type(c) is str

    # test__process_setfl_n_symbols
    eamfile.process_setfl_n_symbols()
    assert type(eamfile.symbols) is list
    assert eamfile.symbols == symbols
    assert type(eamfile.n_symbols) is int
    assert eamfile.n_symbols == n_symbols
    assert eamfile.symbols == ['Ni']
    assert eamfile.n_symbols == 1

    # test__process_setfl_rho_r_args
    eamfile.process_setfl_rho_r_args()

    if False:
        attribute_dict = OrderedDict()
        attribute_dict['symbols'] = eamfile.symbols
        attribute_dict['n_symbols'] = eamfile.n_symbols
        attribute_dict['N_rho'] = eamfile.N_rho
        attribute_dict['d_rho'] = eamfile.d_rho
        attribute_dict['N_r'] = eamfile.N_r
        attribute_dict['d_r'] = eamfile.d_r
        attribute_dict['r_cut'] = eamfile.r_cut
        for k,v in attribute_dict.items():
            print("{}:{}".format(k,v))

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

    #--------------------------------------------------------------------------
    # testing the process_setfl_header_section()
    eamfile = eamtools.EamSetflFile()
    eamfile.read(filename=filename,autoprocess=False)
    eamfile.process_setfl_header_section()

    # this needs to be fixed in the production code
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

    if False:
        attribute_dict = OrderedDict()
        attribute_dict['N_headersection_lines'] = eamfile.N_headersection_lines
        attribute_dict['N_density_lines'] = eamfile.N_density_lines
        attribute_dict['N_embedding_lines'] = eamfile.N_embedding_lines
        attribute_dict['N_pairpotential_lines'] = eamfile.N_pairpotential_lines
        attribute_dict['N_atomicsection_lines'] = eamfile.N_atomicsection_lines
        attribute_dict['N_pairsection_lines'] = eamfile.N_pairsection_lines
        for k,v in attribute_dict.items():
            print('{}:{}'.format(k,v))
        print('n_total_lines:{}'.format(n_total_lines))

    #--------------------------------------------------------------------------
    # testing the process_setfl_atomic_section()
    eamfile = eamtools.EamSetflFile()
    eamfile.read(filename=filename,autoprocess=False)
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

    if False:
        for sym_k,sym_v in eamfile.symbol_dict.items():
            assert type(sym_k) is str
            assert isinstance(sym_v,dict)
            for sym_k_attr,sym_k_value in sym_v.items():
                print("{}:{}:{}".format(
                    sym_k,
                    sym_k_attr,
                    sym_k_value))
        
        for s in eamfile.symbols:
            print('func_embedding[{}]'.format(s))
            print('\ttype:{}'.format(type(eamfile.func_embedding[s])))
            print('\tN:{}'.format(n_row))
            print('func_density[{}]'.format(s))
            print('\ttype:{}'.format(type(eamfile.func_density[s])))
            print('\tN:{}'.format(n_row))
    #--------------------------------------------------------------------------
    # testing the process_setfl_pairpotential_section()
    eamfile = eamtools.EamSetflFile()
    eamfile.read(filename=filename,autoprocess=False)
    eamfile.process_setfl_header_section()
    eamfile.process_setfl_atomic_section()
    eamfile.process_setfl_pairpotential_section()

    assert isinstance(eamfile.func_pairpotential,OrderedDict)
    for pair,pairpot in eamfile.func_pairpotential.items():
        assert isinstance(pairpot,np.ndarray)
        (n_row,n_col) = np.atleast_2d(pairpot).T.shape
        assert eamfile.N_r == n_row


