import numpy as np
import pandas as pd
from pypospack.pyposmat.data import PyposmatDataFile

data_fn = "pyposmat.results.0.out"

data = PyposmatDataFile()
data.read(filename=data_fn)

print(type(data.df))
print(data.df.shape)
print(list(data.df.columns.values))

import matplotlib.pyplot as plt

plt.scatter(
    np.abs(data.df['MgO_NaCl.p11']),
    np.abs(data.df['MgO_NaCl.a0.err'])
)
plt.show()
