
def plot_ellipse(position,covariance,ax=None,**kwargs):
    if ax is None:
        fig, ax = plt.subplots(1,1)

    if covariance.shape == (2,2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1,0],U[0,0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)

    # draw ellipse

    for nsig in range(1,4):
        ax.add_patch(Ellipse(position,nsig*width,nsig*height,angle,**kwargs))

def plot_gmm(gmm,X,label=True,ax=None,dpi=1200,filename=None):
    if ax is None:
        fig, ax = plt.subplots(1,1)
    cluster_id = gmm.fit(X).predict(X)

    if isinstance(X,np.ndarray):
        x = X[:,0]
        y = X[:,1]
    elif isinstance(X,pd.DataFrame):
        x = X[X.columns[0]]
        y = X[X.columns[1]]

    if label:
        ax.scatter(x,y,c=cluster_id,s=1,cmap='viridis',zorder=2)
    else:
        ax.scatter(x,y,s=40,zorder=2)

    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_,gmm.covariances_,gmm.weights_):
        plot_ellipse(pos,covar,alpha=w*w_factor,ax=ax,fill=None)

    ax.set_xlabel(X.columns[0])
    ax.set_ylabel(X.columns[1])
    
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename,dpi=dpi)
