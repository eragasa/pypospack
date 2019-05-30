import pytest
import numpy as np
from z_transformation import ZTransformation

def test____init_____no_args():

    o = ZTransformation()

def dev__fit(X):
    print('type(X)={}'.format(type(X)))
    o = ZTransformation()
    o.fit(X) 
    print('type_checking')
    print('type(o.X)={}'.format(type(o.X)))
    print('type(o.mean)={}'.format(type(o.mean)))
    print('type(o.std)={}'.format(type(o.std)))
    print('X_shape={}'.format(o.X.shape))
    print('mean_shape={}'.format(o.mean.shape))
    print('std_shape={}'.format(o.std.shape))

def dev__fit__w_numpy_array():
    print('dev__fit__w_numpy_array')
    X = np.array([[1,2,3],
                 [2,3,4]])

    dev__fit(X)
def test__fit__w_numpy_array():
    X = np.array([1,2,3])

if __name__ == "__main__":
    dev__fit__w_numpy_array()
