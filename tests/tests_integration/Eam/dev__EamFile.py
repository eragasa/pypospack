import numpy as np
import pypospack.eamtools as eamtools

eamfile = eamtools.EamSetflFile()
comments = ["a","b","c"]
eamfile.comments = comments
for i,v in enumerate(eamfile.comments):
    print('comment[{}]:'.format(str(i)),
          "'{}'".format(eamfile.comments[i]),
          "'{}'".format(comments[i]))

test_setfl_number_format = True
if test_setfl_number_format:
    print(80*"-")
    print("Testing the SETFL number format")

    eamfile = eamtools.EamSetflFile()
    num_format = eamfile.SETFL_NUM_FORMAT
    test_str = 5*num_format.format(1,2,3,45)
    len_test_str = len(test_str)

    print("SETFL_NUM_FORMAT:{}".format(num_format))
    print("str_out:{}".format(test_str))
    print("len_str_out:{}".format(len_test_str))

test_header_section = True

if test_can_read_setfl_file:
    filename = 
if test_header_section:
    print(80*"-")
    print("Testing the header section of the SETFL file")
    
    symbols = ['Ni']
    rho_max = 10
    r_max = 10

    eamfile = eamtools.EamSetflFile()
    eamfile.symbols = symbols
    eamfile.get_comments as string
print(80*"-")
print("Testing the comments section")
eamfile = eamtools.EamSetflFile()
comments = ["comment 1","comment 2","comment 3"]
comment_str = "\n".join(comments)
eamfile.comments = comments
print("TestProduced:\n{}<--".format(eamfile.get_comments_as_string()))
print("Expected:\n{}<--".format(comment_str))
print("Test:",eamfile.get_comments_as_string()==comment_str)

print(80*"-")
print('number of elements line')
eamfile = eamtools.EamSetflFile()
symbols = ['Ni']
tokens = [len(symbols)] + symbols
tokens = [str(s) for s in tokens]
assert isinstance(tokens,list)
n_symbols_line_str = " ".join(tokens)
eamfile.symbols = symbols
print("TestProducted:\n{}<--".format(eamfile.get_n_symbols_line_as_string()))
print("Expected:\n{}<--".format(n_symbols_line_str))
print("Test:",eamfile.get_n_symbols_line_as_string()==n_symbols_line_str)


symbols = ['Ni']
r_max = 10.0
N_r = 500

#r, d_r = eamtools.get_eam_radii_vector(r_max=r_max,N_r=N=r)
#assert isinstance(r,np.ndarray)
#assert isinstance(d_r,float)


