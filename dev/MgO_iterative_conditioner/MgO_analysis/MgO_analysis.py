import numpy as np
import pandas as pd
from pypospack.pyposmat.data import PyposmatDataFile

data_fn = "pyposmat.results.0.out"
#u can turn panda info into numpy array. in matlab this is a pain
#figure out how to do a PCA plot... find examples on scikitlearn python package ( in a different file )
data = PyposmatDataFile()
data.read(filename=data_fn)
#data.df is a pandas dataframe. this is the data structure, How do I get the smallest values? then you can search with this criteria. then I dont
print(type(data.df))
print(data.df.shape)
print(list(data.df.columns.values))
print(data.parameter_names)
print(data.qoi_names)
print(data.error_names)
import matplotlib.pyplot as plt

plt.scatter(
    np.abs(data.df['MgO_NaCl.p11']),
    np.abs(data.df['MgO_NaCl.a0.err'])
)
plt.show()
