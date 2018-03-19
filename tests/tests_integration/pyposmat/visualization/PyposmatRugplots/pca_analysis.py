import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from pypospack.pyposmat.data import PyposmatDataFile

class PcaAnalysis(object):

    def __init__(self):
        self.df = copy.deepcopy(self.df)

        self._configuration = None
        self._df = None

        self.obj_scaler = StandardScaler().fit(self.df)
        self.normed_df = pd.Dataframe(
                self.scaler.transform(self.parameter_df)
                columns=self.parameter_names)
   
    @property
    def parameter_names(self): return self.data.parameter_names

    @property
    def qoi_names(self): return self.data.qoi_names

    @property
    def error_names(self): return self.data.error_names

    @property
    def df(self):
        if self._
    return self._data.df
    
    @df.setter
    def df(self,df):
    @property
    def parameter_df(self): return self.data.df[self.parameter_names]

    @property
    def error_df(self): return self.data.df[self.error_names]

    @property
    def error_df(self): return self.data.df[self.error_names]

    def read_configuration_file(filename):
        self.configuration = PyposmatConfigurationFile
filename = 'subselect.d_metric.bjs.out'

data=PyposmatDataFile(_ilename_bjs)
data.read()

parameter_names = data.parameter_names
qoi_names = data.qoi_names
error_names = data.error_names

parameter_df = data_bjs.df[_parameter_names_bjs]
_qoi_df_bjs = _data_bjs.df[_error_names_bjs]
_error_df_bjs = _data_bjs.df[_error_names_bjs]

_parameter_scaler_bjs = StandardScaler()
_parameter_scaler_bjs.fit(_parameter_df_bjs)
_parameter_norm_bjs = _parameter_scaler_bjs.transform(_parameter_df_bjs)

_parameter_pca_obj_bjs = PCA(n_components='mle',svd_solver='full')
_parameter_pca_bjs = _parameter_pca_obj_bjs.fit_transform(_parameter_norm_bjs)

_param_pca_components = pd.DataFrame(
        _parameter_pca_obj_bjs.components_,
        columns=_parameter_names_bjs)
_param_pca_components['explained_variance'] \
        = _parameter_pca_obj_bjs.explained_variance_
_param_pca_components['explained_variance_ratio'] \
        = _parameter_pca_obj_bjs.explained_variance_ratio_
print(_param_pca_components)

    
nrows, ncols = _parameter_pca_bjs.shape
_parameter_pca_names_bjs = ['parameter_pca_{}'.format(i) for i in range(ncols)]
_parameter_pca_bjs_df = pd.DataFrame(
        data=_parameter_pca_bjs,
        columns=_parameter_pca_names_bjs)

_error_scaler_bjs = StandardScaler()
_error_scaler_bjs.fit(_error_df_bjs)
_error_norm_bjs = _error_scaler_bjs.transform(_error_df_bjs)

_error_pca_obj_bjs = PCA()
_error_pca_bjs = _error_pca_obj_bjs.fit_transform(_error_norm_bjs)
nrows, ncols = _error_pca_bjs.shape
_error_pca_names_bjs = ['error_pca_{}'.format(i) for i in range(ncols)]
_error_pca_bjs_df = pd.DataFrame(
        data=_error_pca_bjs,
        columns=_error_pca_names_bjs)
print(_parameter_pca_names_bjs)
print(_error_pca_names_bjs)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig_bjs = plt.figure(figsize=plt.figaspect(0.5))
ax_param_bjs = fig_bjs.add_subplot(1,2,1,projection='3d')
ax_param_bjs.scatter(
        _parameter_pca_bjs_df['parameter_pca_1'],
        _parameter_pca_bjs_df['parameter_pca_2'],
        _parameter_pca_bjs_df['parameter_pca_3'],
        marker='.',
        s=1)
ax_error_bjs = fig_bjs.add_subplot(1,2,2,projection='3d')
ax_error_bjs.scatter(
        _error_pca_bjs_df['error_pca_1'],
        _error_pca_bjs_df['error_pca_2'],
        _error_pca_bjs_df['error_pca_3'],
        marker='.',
        s=1)
fig_bjs.savefig('pca_bjs.png')

_data_fs=PyposmatDataFile(_filename_fs)
_data_fs.read()

_parameter_names_fs = _data_fs.parameter_names
_qoi_names_fs = _data_fs.qoi_names
_error_names_fs = _data_fs.error_names

_parameter_df_fs = _data_fs.df[_parameter_names_fs]
_qoi_df_fs = _data_fs.df[_error_names_fs]
_error_df_fs = _data_fs.df[_error_names_fs]

_parameter_scaler_fs = StandardScaler()
_parameter_scaler_fs.fit(_parameter_df_fs)
_parameter_norm_fs = _parameter_scaler_fs.transform(_parameter_df_fs)

_parameter_pca_obj_fs = PCA()
_parameter_pca_fs = _parameter_pca_obj_fs.fit_transform(_parameter_norm_fs)
nrows, ncols = _parameter_pca_fs.shape
_parameter_pca_names_fs = ['parameter_pca_{}'.format(i) for i in range(ncols)]
_parameter_pca_fs_df = pd.DataFrame(
        data=_parameter_pca_fs,
        columns=_parameter_pca_names_fs)

_error_scaler_fs = StandardScaler()
_error_scaler_fs.fit(_error_df_fs)
_error_norm_fs = _error_scaler_fs.transform(_error_df_fs)

_error_pca_obj_fs = PCA()
_error_pca_fs = _error_pca_obj_fs.fit_transform(_error_norm_fs)
nrows, ncols = _error_pca_fs.shape
_error_pca_names_fs = ['error_pca_{}'.format(i) for i in range(ncols)]
_error_pca_fs_df = pd.DataFrame(
        data=_error_pca_fs,
        columns=_error_pca_names_fs)
print(_parameter_pca_names_fs)
print(_error_pca_names_fs)

fig_fs = plt.figure(figsize=plt.figaspect(0.5))
ax_param_fs = fig_fs.add_subplot(1,2,1,projection='3d')
ax_param_fs.scatter(
        _parameter_pca_fs_df['parameter_pca_1'],
        _parameter_pca_fs_df['parameter_pca_2'],
        _parameter_pca_fs_df['parameter_pca_3'],
        marker='.',
        s=1)
ax_error_fs = fig_fs.add_subplot(1,2,2,projection='3d')
ax_error_fs.scatter(
        _error_pca_fs_df['error_pca_1'],
        _error_pca_fs_df['error_pca_2'],
        _error_pca_fs_df['error_pca_3'],
        marker='.',
        s=1)
fig_fs.savefig('pca_fs.png')

_error_df = pd.concat([_error_df_bjs,_error_df_fs])
_error_scaler = StandardScaler()
_error_scaler.fit(_error_df)
_error_norm_bjs = _error_scaler.transform(_error_df_bjs)
_error_norm_fs = _error_scaler.transform(_error_df_fs)

_error_pca_obj = PCA()
_error_pca_fs = _error_pca_obj.fit_transform(_error_norm_fs)
_error_pca_bjs = _error_pca_obj.fit_transform(_error_norm_bjs)
nrows, ncols = _error_pca_fs.shape
_error_pca_names = ['error_pca_{}'.format(i) for i in range(ncols)]
_error_pca_fs_df = pd.DataFrame(
        data=_error_pca_fs,
        columns=_error_pca_names)
_error_pca_bjs_df = pd.DataFrame(
        data=_error_pca_bjs,
        columns=_error_pca_names)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax_error = fig_fs.add_subplot(1,1,1,projection='3d')
ax_error.scatter(
        _error_pca_fs_df['error_pca_1'],
        _error_pca_fs_df['error_pca_2'],
        _error_pca_fs_df['error_pca_3'],
        marker='.',
        s=1,
        c='b')
ax_error.scatter(
        _error_pca_bjs_df['error_pca_1'],
        _error_pca_bjs_df['error_pca_2'],
        _error_pca_bjs_df['error_pca_3'],
        marker='.',
        s=1,
        c='r')
fig_fs.savefig('pca.png')

#http://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html#sphx-glr-auto-examples-preprocessing-plot-all-scaling-py
