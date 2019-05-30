import numpy as np

class ZTransformation(object):

    def __init__(self,X=None,mean=None,std=None):

        self.mean = None
        self.std = None

    def fit(self,X):
        self.X = X
        self.mean = X.mean()
        self.std =X.std()

    def transform(self,X):
        pass  
if __name__ == "__main__":

    a = [1,2,3]
    a_mean = np.array(a).mean()
    print(a_mean
