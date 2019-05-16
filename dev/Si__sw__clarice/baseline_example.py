"""
tutorial file on basic file manipulation
"""
import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

# this should work if your paths are correct
# we have established this method of collaboration so that your scripts will work across
# platforms since all data is in the expected directory, and doesn't have to changed 
# for each individuals
pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
print("pypospack_root_dir={}".format(pypospack_root_dir))

# path to the configuration file
config_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.config.in')
print("config_fn={}".format(config_fn))

# path to the data file
data_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')
print("data_fn={}".format(data_fn))

#read the configuration file into an instance object
o_config = PyposmatConfigurationFile()
o_config.read(filename=config_fn)

#parameter_names
print('parameter_names')
print(type(o_config.parameter_names))
for k in o_config.parameter_names:
    print('\t',k)

#qoi_names
print('qoi_names')
print(type(o_config.qoi_names))
for k in o_config.qoi_names:
    print('\t',k)

#read the datafile into an instance of an object
o_data = PyposmatDataFile()
o_data.read(filename=data_fn)

# df is a pandas dataframe attribute
print(type(o_data.df))

# selecting only the qoi columns
print(o_data.df[o_config.parameter_names])

# selecting only the parameter columns
print(o_data.df[o_config.qoi_names])
