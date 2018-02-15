import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from real_data_pca import read_data_file

def compare_inversion(df):
    # start with raw data df and do transform and inversion in here
    # remove sim_id column for analysis
    orig = df.drop(['sim_id'], axis=1)
    # write raw df to csv
    orig.to_csv(path_or_buf='raw_data.csv', sep=',')
    orig_np, orig_norms = normalize(orig, return_norm=True)
    nrows, ncols = orig_np.shape
    orig_names = ['orig_{}'.format(i) for i in range(ncols)]
    orig_df = pd.DataFrame(data=orig_np, columns=orig_names)

    obj_pca = PCA()
    pca_np = obj_pca.fit_transform(orig_df)
    pca_names = ['pca_{}'.format(i) for i in range(ncols)]
    pca_df = pd.DataFrame(data=pca_np, columns=pca_names)

    # undo pca to get to normalized data
    inv_np = obj_pca.inverse_transform(pca_np)
    inv_names = ['inv_{}'.format(i) for i in range(ncols)]
    inv_df = pd.DataFrame(data=inv_np, columns=inv_names)

    assert orig_df.shape == pca_df.shape == inv_df.shape

    # undo normalization to get to raw data
    for i, n in enumerate(orig_norms.tolist()):
        inv_df.iloc[i] = inv_df.iloc[i] * n

    # write denormalized df to csv
    inv_df.to_csv(path_or_buf='fully_inverted_data.csv', sep=',')

    # compare error between raw input and inverted output
    difference = np.array(inv_df) - np.array(orig)

    norm = stats.norm(scale=difference.std(), loc=difference.mean())
    x = np.linspace(difference.min(), difference.max(), 100)
    plt.plot(x, norm.pdf(x), 'r-', lw=5, alpha=0.6, label='pdf of errors')
    plt.hist(difference)
    plt.legend(loc='best')
    plt.title(s='Error Analysis of PCA Inversion')
    plt.xlabel(s='error')
    plt.ylabel(s='frequency')
    plt.show()


if __name__ == "__main__":
    results = read_data_file(r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out')
    compare_inversion(df=results[0]['total_df'])
