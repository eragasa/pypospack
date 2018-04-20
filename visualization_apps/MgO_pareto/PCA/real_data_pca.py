import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering, AffinityPropagation
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler


def read_data_file(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    names = [s.strip() for s in lines[0].strip().split(',')]
    types = [s.strip() for s in lines[1].strip().split(',')]

    all_values = []
    for i in range(2, len(lines)):
        line = lines[i].strip()
        values = [float(s.strip()) for s in line.split(',')]
        values[0] = int(values[0])
        all_values.append(list(values))
    values = np.array(all_values)

    parameter_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'param']
    qoi_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'qoi']
    error_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'err']
    df = pd.DataFrame(data=values, columns=names, copy=True)
    df.set_index('sim_id')
    parameter_df = df[parameter_names]
    # static params are: chrg_O, MgMg_A, MgMg_rho, MgMg_C, MgO_C
    parameter_df = parameter_df.drop(['chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_C'], axis=1)
    # make rho exponential to put on equal standing with other factors
    parameter_df['MgO_rho'] = np.exp(parameter_df['MgO_rho'])
    parameter_df['OO_rho'] = np.exp(parameter_df['OO_rho'])
    error_df = df[error_names]
    qoi_df = df[qoi_names]

    return {'total_df': df, 'param_df': parameter_df, 'err_df': error_df, 'qoi_df': qoi_df}, \
           {'param_names': parameter_names, 'err_names': error_names, 'qoi_names': qoi_names}


class TransformPCA(object):

    def __init__(self, num_clusters, param_df, err_df, qoi_df, total_df, cluster_by='kmeans'):
        self.num_clusters = num_clusters
        self.cluster_by = cluster_by
        '''
        change input method to only a filepath:
        with open(filepath) as f:
            # data file class must be modifed to accept this
            data = PyposmatDataFile(f)
        self._param_df = data.parameter_df
        self._err_df = data.error_df
        self._qoi_df = data.qoi_df
        self._total_df = data.df
        assert self._total_df is not None
        '''
        self.param_df = param_df
        self.err_df = err_df
        self.qoi_df = qoi_df
        self.total_df = total_df
        self.pca_df = pd.DataFrame()
        self._normalize_and_transform()

    def _normalize_and_transform(self):
        pca_df_collection = []
        clusters = None
        for df, label in zip([self.param_df, self.err_df, self.qoi_df], ['param', 'err', 'qoi']):
            # temporary abs value test
            if label == 'err':
                df = np.absolute(df)
            # ------------------------
            df = normalize(df)
            #norm = StandardScaler()
            #df = norm.fit_transform(df)
            obj_pca = PCA()
            pca_np = obj_pca.fit_transform(df)
            nrows, ncols = pca_np.shape
            pca_names = [label + '_pca_{}'.format(i) for i in range(ncols)]
            pca_df = pd.DataFrame(data=pca_np, columns=pca_names)
            pca_df_collection.append(pca_df)
            self._write_excel(obj_pca, pca_names, label)
            # use the param data to define clusters
            if label == 'param':
                if self.cluster_by == 'dbscan':
                    clusters = self._dbscan_cluster(pca_df)
                elif self.cluster_by == 'agg':
                    clusters = self._agglomerative_cluster(pca_df)
                elif self.cluster_by == 'aff_pro':
                    clusters = self._affinity_cluster(pca_df)
                else:
                    clusters = self._kmeans_cluster(pca_df)
        self.pca_df = pd.DataFrame(pd.concat(pca_df_collection, axis=1))
        self.pca_df['cluster_id'] = clusters

    def _write_excel(self, obj_pca, pca_names, label):
        # !!! The labels this function produces are wrong !!!
        # be sure to correct by hand for now
        eigen_values = obj_pca.explained_variance_
        eigen_vectors = obj_pca.components_
        eigen_df = pd.DataFrame(columns=pca_names)
        eigen_vectors = np.transpose(eigen_vectors)
        # .tolist() actually returns iterable class not list class
        eigen_vectors = list(eigen_vectors.tolist())
        eigen_df['eigen_values'] = eigen_values
        for i, names in enumerate(pca_names):
            eigen_df[names] = eigen_vectors[i]
        writer = pd.ExcelWriter(label + '_eigen_table.xlsx')
        eigen_df.to_excel(writer, label)
        writer.save()

    def _kmeans_cluster(self, param_pca_df):
        # do pca on only param and get cluster groupings for that
        obj_kmeans = KMeans(n_clusters=self.num_clusters)
        obj_kmeans = obj_kmeans.fit(param_pca_df)
        labels = obj_kmeans.labels_
        return labels

    def _dbscan_cluster(self, param_pca_df):
        obj_dbscan = DBSCAN(eps=0.3, min_samples=100)
        obj_dbscan = obj_dbscan.fit(param_pca_df)
        labels = obj_dbscan.labels_
        return labels

    def _agglomerative_cluster(self, param_pca_df):
        obj_agg = AgglomerativeClustering(n_clusters=self.num_clusters)
        obj_agg = obj_agg.fit(param_pca_df)
        labels = obj_agg.labels_
        return labels

    def _affinity_cluster(self, param_pca_df):
        obj_aff = AffinityPropagation(damping=0.8, convergence_iter=10)
        obj_aff = obj_aff.fit(param_pca_df)
        labels = obj_aff.labels_
        return labels

# ^^^-------- Private methods to do stuff
# vvv-------- Public methods to visualize

    def plot_kde_2d(self, dtype):
        fig = plt.figure(figsize=(8,6))
        _clusterids = set(self.pca_df['cluster_id'])
        for clusterid in _clusterids:
            vect0 = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_0']
            vect1 = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_1']
            xmin, xmax = vect0.min(), vect0.max()
            ymin, ymax = vect1.min(), vect1.max()

            # Peform the kernel density estimate
            xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
            positions = np.vstack([xx.ravel(), yy.ravel()])
            values = np.vstack([vect0, vect1])
            kernel = stats.gaussian_kde(values)
            f = np.reshape(kernel(positions).T, xx.shape)

            ax = fig.add_subplot(311+clusterid)
            # Contourf plot
            ax.contourf(xx, yy, f, cmap='Blues')
            # Contour plot
            cset = ax.contour(xx, yy, f, colors='k')
            # Label plot
            ax.clabel(cset, inline=1, fontsize=10)
            ax.set_xlabel('PCA 0')
            ax.set_ylabel('PCA 1')
            ax.set_title(label=dtype+' Group '+str(clusterid)+' Distribution')
        plt.subplots_adjust(hspace=1.0)
        plt.show()

    def plot_2d(self, dtype):
        _clusterids = set(self.pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for clusterid in _clusterids:
            x = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype + '_pca_0']
            y = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype + '_pca_1']
            ax.scatter(x, y, s=1)
        plt.title(s=dtype + ' PCA Transformation')
        plt.show()

    def plot_3d(self, dtype):
        _clusterids = set(self.pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for clusterid in _clusterids:
            x = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_0']
            y = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_1']
            z = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_2']
            print(x.shape, y.shape, z.shape)
            ax.scatter(x, y, z, s=1)
        plt.title(s=dtype+' PCA Transformation')
        plt.show()

    # replaced by the kde plot
    def plot_hist(self, dtype):
        fig = plt.figure(figsize=(8,6))
        _clusterids = set(self.pca_df['cluster_id'])
        # make each axis a different subplot
        # each group is represented on each subplot
        for dim in range(3):    # pca0, pca1, pca2
            ax = fig.add_subplot(310+(dim+1))
            for clusterid in _clusterids:   # group1, group2, group3
                vect = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][dtype+'_pca_'+str(dim)]
                ax.hist(vect, bins=25)
                ax.set_title("PCA Axis "+str(dim), x=0)
        fig.tight_layout(h_pad=1)
        fig.suptitle(t=dtype+" Cluster Distribution")
        plt.show()

if __name__ == "__main__":
    # fix filepath issue
    '''
    sys.path.append('../visualization_apps/MgO_pareto/')
    data_folder = 'data'
    filename = 'culled_009.out'
    data = PyposmatDataFile(filename=os.path.join(data_folder, filename))
    _param_df = data.parameter_df
    _err_df = data.error_df
    _qoi_df = data.qoi_df
    _total_df = data.df
    tp = TransformPCA(num_clusters=3, param_df=_param_df, err_df=_err_df, qoi_df=_qoi_df, total_df=_total_df)
    '''

    # cluster_by kmeans, dbscan, agg, or aff_pro (kmeans default)
    results = read_data_file(r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out')
    tp = TransformPCA(num_clusters=3, param_df=results[0]['param_df'],
                      err_df=results[0]['err_df'], qoi_df=results[0]['qoi_df'],
                      total_df=results[0]['total_df'],
                      cluster_by='kmeans')
    # dtype must be 'param', 'err', or 'qoi'
    tp.plot_3d(dtype='param')
    tp.plot_3d(dtype='err')
