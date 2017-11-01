import pypospack.eamtools as eamtools
import pypospack.io.filesystem as filesystem
from collections import OrderedDict
import numpy as np

class DevEamSetflFile(eamtools.EamSetflFile):
    
    def write(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.filename = filename

    def read(self,filename=None,autoprocess=True):
        if filename is not None:
            assert type(filename) is str
            self.filename = filename

        self.lines = filesystem.read_file_as_lines(self.filename)

        assert isinstance(self.lines,list)
        if autoprocess is True:
            self.process_setfl_header_section()


if __name__ == "__main__":
    import os

    filename = '../Ni1_Mendelev_2010.eam.fs.txt'
    assert os.path.isfile(filename)

    eamfile = DevEamSetflFile()
    print("N_headersection_lines:{}".format(
        eamfile.N_headersection_lines))
    eamfile.read(filename=filename)

    print("eam_filename:{}".format(
        filename))
    print("#"+79*"-")
    print("{:^79}".format(
        "comment section"))
    print("#"+79*"-")
    print(eamfile.get_comments_as_string())
    assert isinstance(eamfile.comments,list)
    assert len(eamfile.comments) == 3
    
    print("#"+79*"-")
    print("{:^79}".format(
        "some internal attributes"))
    print("#"+79*"-")
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

    print("#"+79*"-")
    print("{:^79}".format(
        "symbol information"))
    print("#"+79*"-")
    eamfile = DevEamSetflFile()
    eamfile.read(filename=filename,autoprocess=False)
    eamfile.process_setfl_header_section()
    eamfile.process_setfl_atomic_section()
    eamfile.process_setfl_pairpotential_section()
    for sym_k,sym_v in eamfile.symbol_dict.items():
        assert type(sym_k) is str
        assert isinstance(sym_v,dict)
        for sym_k_attr,sym_k_value in sym_v.items():
            print("{}:{}:{}".format(
                sym_k,
                sym_k_attr,
                sym_k_value))

    for s in eamfile.symbols:
        (n_row,n_col) = np.atleast_2d(eamfile.func_embedding[s]).T.shape
        print('func_embedding[{}]'.format(s))
        print('\ttype:{}'.format(type(eamfile.func_embedding[s])))
        print('\tN:{}'.format(n_row))
        (n_row,n_col) = np.atleast_2d(eamfile.func_embedding[s]).T.shape
        print('func_density[{}]'.format(s))
        print('\ttype:{}'.format(type(eamfile.func_density[s])))
        print('\tN:{}'.format(n_row))
