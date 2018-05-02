import numpy as np
import pandas as pd
from pypospack.pyposmat.data import PyposmatDataFile

data_fn = "pyposmat.results.0.out"
#u can turn panda info into numpy array. in matlab this is a pain
#figure out how to do a PCA plot... find examples on scikitlearn python package
#( in a different file )
data = PyposmatDataFile()
data.read(filename=data_fn)
#look up how to do kernel density estimate.. heat maps are cool!!!!!!!!
#wiki.materialsexmachina.com/index.php/Kernel_Density_Estimate
#youre chosing axes that a re linear ocmbinations of the parameters.
#youre also choosing the first axes/second axes that are orthogonal to each
#other and at the first axis is the direction of the greatest variace and the
#second axes is in the orthogonal direction which describes  the greatest amount
# of variance in a direction orthogonal to the first PCA vector.
#data.df is a pandas dataframe. this is the data structure, How do I get the
#smallest values? then you can search with this criteria. then I dont

#Data information

print('Data structure =')
print(type(data.df))
print('Shape =')
print(data.df.shape)
print('Columns = ')
print(list(data.df.columns.values))
print('Parameter names =')
print(data.parameter_names)
print('QOI names =')
print(data.qoi_names)
print('Error names =')
print(data.error_names)

#Plotting
import matplotlib.pyplot as plt

#print(data.df.shape)
#df = data.df.loc[(
#df['MgO_NaCl.p11'] < 1000) &
#df['MgO_NaCl.p11'] > 0) &
#df['MgO_NaCl.a0.err'] > 0) &
#df['MgO_NaCl.a0.err'] < 1000)
#]
#print(df.shape)
# df_xlim=data.df['MgO_NaCl.p11']<10000
# for x in
# if df_xlim == True
# truncated_x=[truncated_x ]
#
# print(df_xlim)
# print(df_xlim.shape)
plt.scatter(
    np.abs(data.df['MgO_NaCl.p11']/1000),
    np.abs(data.df['MgO_NaCl.a0.err'])
)
# plt.xlim(-10000,  10000)
# plt.ylim(-10000, 10000)
plt.title("Error of Pressures on MgO Cubic Cell")
plt.xlabel("NaCl p11 Tensor (kilobar)")
plt.ylabel("Error")
plt.show()
