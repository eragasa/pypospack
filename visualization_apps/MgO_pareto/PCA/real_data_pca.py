import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize
from sklearn.neighbors.kde import KernelDensity

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
    error_df = df[error_names]
    qoi_df = df[qoi_names]

    return {'total_df': df, 'param_df': parameter_df, 'err_df': error_df, 'qoi_df': qoi_df}, \
           {'param_names': parameter_names, 'err_names': error_names, 'qoi_names': qoi_names}


class TransformPCA(object):

    def __init__(self, num_clusters, param_df, err_df, qoi_df, total_df):
        self.num_clusters = num_clusters
        self.param_df = param_df
        self.err_df = err_df
        self.qoi_df = qoi_df
        self.total_df = total_df
        self.pca_df = pd.DataFrame()
        self.normalize_and_transform()

    def normalize_and_transform(self):
        pca_df_collection = []
        clusters = None
        for df, label in zip([self.param_df, self.err_df, self.qoi_df], ['param', 'err', 'qoi']):
            df = normalize(df)
            obj_pca = PCA()
            pca_np = obj_pca.fit_transform(df)
            nrows, ncols = pca_np.shape
            pca_names = [label + '_pca_{}'.format(i) for i in range(ncols)]
            pca_df = pd.DataFrame(data=pca_np, columns=pca_names)
            pca_df_collection.append(pca_df)
            # self.write_excel(obj_pca, pca_names, label)
            # use the param data to define clusters
            if label == 'param':
                clusters = self.kmeans_cluster(pca_df)
        self.pca_df = pd.DataFrame(pd.concat(pca_df_collection, axis=1))
        self.pca_df['cluster_id'] = clusters

    def write_excel(self, obj_pca, pca_names, label):
        # create df of vectors and values to write excel file easily
        eigen_values = obj_pca.explained_variance_
        eigen_vectors = obj_pca.components_
        eigen_df = pd.DataFrame(columns=pca_names)
        eigen_vectors = eigen_vectors.tolist()
        for i, names in enumerate(pca_names):
            eigen_df[names] = eigen_vectors[i]
            if i == len(eigen_vectors) - 1:
                eigen_df.loc[i + 1] = eigen_values.tolist()
        writer = pd.ExcelWriter(label + '_eigen_table.xlsx')
        eigen_df.to_excel(writer, label)
        writer.save()

    def kmeans_cluster(self, param_pca_df):
        # do pca on only param and get cluster groupings for that
        obj_kmeans = KMeans(n_clusters=self.num_clusters)
        obj_kmeans = obj_kmeans.fit(param_pca_df)
        labels = obj_kmeans.labels_
        return labels

    def plot_2d(self):
        _clusterids = set(self.pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for label in ['param', 'err', 'qoi']:
            for clusterid in _clusterids:
                x = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][label + '_pca_0']
                y = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid][label + '_pca_1']
                ax.scatter(x, y, s=1)
            plt.title(s=label + ' PCA transformation')
            plt.show()

    def plot_3d(self):
        _clusterids = set(self.pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for clusterid in _clusterids:
            x = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['param_pca_0']
            y = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['param_pca_1']
            z = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['param_pca_2']
            print(x.shape, y.shape, z.shape)
            ax.scatter(x, y, z, s=1)
        plt.title(s='Parameter PCA')
        plt.show()

    def plot_hist(self):
        fig = plt.figure(figsize=(8,6))
        _clusterids = set(self.pca_df['cluster_id'])
        # make each axis a different subplot
        # each group is represented on each subplot
        for dim in range(3):    # pca0, pca1, pca2
            ax = fig.add_subplot(int('31'+str(dim+1)))
            for clusterid in _clusterids:   # group1, group2, group3
                vect = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['param_pca_'+str(dim)]
                ax.hist(vect, bins=25)
                ax.set_title("PCA Axis "+str(dim), x=0)
        fig.tight_layout(h_pad=1)
        fig.suptitle(t="Param Cluster Distribution")
        plt.show()

if __name__ == "__main__":
    results = read_data_file(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')
    tp = TransformPCA(num_clusters=3, param_df=results[0]['param_df'],
                      err_df=results[0]['err_df'], qoi_df=results[0]['qoi_df'],
                      total_df=results[0]['total_df'])
    tp.plot_hist()
