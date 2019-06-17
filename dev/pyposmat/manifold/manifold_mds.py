from manifold import Manifold

class MdsManifold(Manifold):

    manifold_type = 'MDS'

    def __init__(self,pyposmat_configuration,pyposmat_data,manifold_config=None):
        Manifold.__init__(self,
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat_data,
                manifold_config=manifold_config
                )

    def initialize_manifold_configuration(self,manifold_configuration=None):
        DEFAULT_MDS_MANIFOLD_CONFIGURATION = OrderedDict()
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['n_components'] = 2
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['init'] = 'pca'
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['random_state'] = 0

        if manifold_configuration is None:
            self.manifold_configuration = DEFAULT_TSNE_MANIFOLD_CONFIGURATION
        elif isinstance(manifold_configuration,dict):
            self.manifold_configuration = manifold_configuration

            for k,v in DEFAULT_MDS_MANIFOLD_CONFIGURATION.items():
                if k not in self.manifold_configuration:
                    self.manifold_configuration[k] = v

    def learn_manifold(self,names,scaling_type='standard'):
        if isinstance(names,list):
            names_ = names
        elif names == 'qois':
            names_ = self.configuration.qoi_names
        elif names == 'free_parameters':
            names_ = self.configuration.free_parmaeter_names
        elif names == 'all':
            names_ = self.configuration.qoi_names \
                    + self.configuration.free_parameter_names

        X = self.data.df[names_]
        X_scaled = get_scaled_data(X,scaling_type=scaling_type)
        self.manifold = manifold.MDS(**self.manifold_configuration)
        Y = mds.fit_transform(X_scaled)
        
        for i in self.manifold_configuration['n_components']:
            self.data.df['MDS_{}'.format(str(i))] = Y[:,i]

    def MDS(self,i):
        return self.data.df['MDS_{}'.format(str(i))]

    def scale_data(self,X,scaling_type):
        if scaling_type='standard':
            X_scaled = preprocessing.scale(X)
        elif scaling_type='none':
            X_scaled = X
        else:
            raise ValueError('unknown scaling type')
