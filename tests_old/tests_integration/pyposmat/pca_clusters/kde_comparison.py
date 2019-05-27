import os
from collections import OrderedDict
from scipy.stats import gaussian_kde
from pypospack.statistics import kullbach_lieber_divergence
from pypospack.kde import Chiu1999_h,Silverman1986_h 
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile



src_dir = "../../../data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_02"
config_fn = os.path.join(src_dir,'pyposmat.config.in')

config = PyposmatConfigurationFile()
config.read(filename=config_fn)
n_iterations = config.n_iterations

# load data
data = OrderedDict()
data['kde'] = OrderedDict()
print('loading datafiles')
for i in range(n_iterations):
    data_fn = os.path.join(src_dir,'pyposmat.kde.{}.out'.format(i))
    print('reading {}...'.format(data_fn))
    data['kde'][i] = PyposmatDataFile()
    try:
        data['kde'][i].read(filename=data_fn)
    except FileNotFoundError as e:
        print("the number of max iterations is actually {}".format(i-1))
        n_iterations = i-1
        break

# comparing kde estimates
print('compare kde estimates')
kld = [-1]
for i in range(1,n_iterations):
    df_0 = data['kde'][i-1].df
    df_1 = data['kde'][i].df
    df_0_p = df_0[config.parameter_names]
    df_1_p = df_1[config.parameter_names]
    nr0,nc0= df_0_p.shape
    nr1,nc1= df_1_p.shape
    #print('nrows:: {}={},{}={}'.format(i,nr0,i+1,nr1))
    #print('ncols:: {}={},{}={}'.format(i,nc0,i+1,nc1))
    X0 = df_0_p.values
    X1 = df_1_p.values
    #print('X0:',X0.shape,type(X0))
    silverman86_h0 = Silverman1986_h(X0.T)
    silverman86_h1 = Silverman1986_h(X1.T)
    chiu99_h0 = Chiu1999_h(X0.T)
    chiu99_h1 = Chiu1999_h(X1.T)
    kde_0=gaussian_kde(X0.T,chiu99_h0)
    kde_1=gaussian_kde(X1.T,chiu99_h1)
    kld.append(kullbach_lieber_divergence(kde_0,kde_1,400))
    print(i,silverman86_h1,chiu99_h1,kld[i])
for i,v in enumerate(kld):
    print(i,v)
