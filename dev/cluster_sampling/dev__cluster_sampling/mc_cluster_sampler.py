
class PyposmatClusterSampler(PyposmatEngine):

    def __init__(self,
            filename_in='pyposmat.config.in',
            filename_out='pyposmat.results.out',
            mpi_rank=None,
            mpi_size=None,
            base_directory=None):

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)
        self.mpi_rank=mpi_rank
        self.mpi_size=mpi_size
        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out
        self.pyposmat_data_bad_filename = 'pypospack.results.bad'
        # for cluster analysis
        self.pyposmat_configuration_fn = filename_in

    @property
    def configuration_fn(self): return self.pyposmat_configuration_fn

    @property
    def data_fn(self): return self.pyposmat_data_in_filename

    def configure_pyposmat_datafile_in(self,filename):
        self.pyposmat_data_in_filename = filename
        self.pyposmat_datafile_in = PyposmatDataFile(filename)
        self.data = PyposmatDataFile()
        self.data.read(filename)

    def configure_pyposmat_datafile_out(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.pyposmat_data_out_filename = filename
        self.pyposmat_datafile_out = PyposmatDataFile(filename)

    def read_configuration_file(self,filename=None):
        if filename is not None:
            _filename = filename
            self.pyposmat_configuration_fn = filename
        else:
            _filename = self.pyposmat_configuration_fn

        try:
            PyposmatEngine.read_configuration_file(self,filename=_filename)
        except FileNotFoundError as e:
            print("Cannot read filename:")
            print("    filename={}".format(filename))
            print("    o.pyposmat_configuration_fn={}".format(self.pyposmat_configuration_fn))
            print("    _filename={}".format(_filename))
            raise e

        # get information from the PyposmatConfigurationFile object
        self.structure_directory = self.configuration.structures['structure_directory']
        self.parameter_names = [p for p in self.configuration.sampling_distribution]
        self.qoi_names = [k for k in self.configuration.qois]
        self.error_names = ['{}.err'.format(k) for k in self.qoi_names]
        self.parameter_distribution_definition =\
                self.configuration.sampling_distribution
        
        try:
            self.free_parameter_names = [k for k,v in self.parameter_distribution_definition.items() if v[0] != 'equals']
        except KeyError as e:
            print(self.parameter_distribution_definition.items())
            raise

        if self.configuration.sampling_constraints is not None:
            self.parameter_constraints = copy.deepcopy(self.configuration.sampling_constraints)
        else:
            self.parameter_constraints = OrderedDict()

        self.constrained_parameter_names = []
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                self.constrained_parameter_names.append(p)

    def run_simulations(self,
            i_iteration,
            n_samples=None,
            filename=None):
        # process the arguments of the method first

        # arg: i_iteration
        i = i_iteration

        # arg: n_samples
        if n_samples is not None:
            _n_samples = n_samples
        else:
            try:
               _n_samples = self.configuration.sampling_type[i]['n_samples_per_cluster']
            except KeyError as e:
                print("must use \"n_samples_per_cluster\" keyword to describe the number of simulations per cluster")
                raise e

        # arg: filename
        if filename is None:
            _filename = self.pyposmat_data_in_filename
        else:
            _filename = filename

        # read in the the datafile
        self.data = PyposmatDataFile()
        self.data.read(_filename)

        # determine the sampling type
        _sampling_type = self.configuration.sampling_type[i]['type']
        if _sampling_type == 'kde_w_clusters':

            if 'cluster_id' in self.data.df.columns:
                print("i still don't do anything, and it's eugene's fault")
                pass
                # Eugene needs to write some code for this to work in PyposmatDataFile and ClusterAnalysis
            else:
                print('clustering data from:\n{}'.format(_filename))
                # prepare the clustering params in an ordered dict
                
                cluster_args = None
                if 'cluster_args' not in self.configuration.sampling_type[i]:
                    cluster_args = self.get_clustering_parameters(
                            configuration_fn=self.configuration_fn,
                            data_fn=_filename)
                else:
                    cluster_args = self.configuration.sampling_type[i]['cluster_args']

                # send the data to be clustered
                obj_cluster_analysis = PyposmatClusterAnalysis.init_from_ordered_dict(cluster_args)
                obj_cluster_analysis.preprocess_data(cluster_args)
                obj_cluster_analysis.calculate_manifold(cluster_args)
                obj_cluster_analysis.calculate_kNN_analysis(cluster_args)
                obj_cluster_analysis.calculate_clusters(cluster_args)
                
                # use newly clustered data for sampling
                self.data.df = obj_cluster_analysis.data.df
                
                print(self.data.df)
                print("if i do or do not do anything it's seaton's fault")

        else:
            raise ValueError(
                'unknown sampling type:{}'.format(
                    _sampling_type
                )
            )
        # actual usage of the clustered data goes here
        # get unique cluster ids
        cluster_ids = set(self.data.df['cluster_id'])
        for cluster_id in cluster_ids:
            print(cluster_id)
            mc_sampler = PyposmatMonteCarloSampler()
            mc_sampler.pyposmat_filename_out = 'iammadeofmagic.out'
            mc_sampler.configuration = PyposmatConfigurationFile()
            mc_sampler.configuration.read(self.configuration_fn)
            mc_sampler.configuration.sampling_type[i] = OrderedDict()
            mc_sampler.configuration.sampling_type[i]['type'] = 'kde'
            mc_sampler.configuration.sampling_type[i]['n_samples'] = _n_samples
            
            mc_sampler.create_base_directories()
            mc_sampler.read_configuration_file(self.configuration_fn)
            mc_sampler.configure_qoi_manager()
            mc_sampler.configure_task_manager()
            mc_sampler.configure_pyposmat_datafile_out()
            #pyposmat_datafile_out = PyposmatDataFile(filename_out)


            mc_sampler.print_structure_database()
            mc_sampler.print_sampling_configuration()
            mc_sampler.print_initial_parameter_distribution()

            for i in range(mc_sampler.n_iterations):
                #engine.run_simulations(i_iteration=i,n_samples=1000)
                mc_sampler.run_kde_sampling(
                    n_samples=_n_samples,
                    filename_in=self.data_fn,
                    cluster_id=cluster_id
                )  

        if self.mpi_rank == 0:
            print(i_iteration,_n_samples,_sampling_type)

    def get_clustering_parameters(
            self,
            configuration_fn,
            data_fn
            ):

        d = OrderedDict()
        d['configuration_fn'] = configuration_fn
        d['data_fn'] = data_fn
        d['include_parameters'] = True
        d['include_qois'] = True
        d['include_errors'] = False

        d['preprocessing'] = OrderedDict()
        d['preprocessing']['type'] = 'standard_scaler'
        d['preprocessing']['args'] = OrderedDict()
        d['preprocessing']['args']['copy'] = True
        d['preprocessing']['args']['with_mean'] = True
        d['preprocessing']['args']['with_std'] = True

        d['manifold'] = OrderedDict()
        d['manifold']['type'] = 'tsne'
        d['manifold']['args'] = OrderedDict()
        d['manifold']['args']['n_components'] = 2
        d['manifold']['args']['perplexity'] = 30
        d['manifold']['args']['early_exaggeration'] = 12
        d['manifold']['args']['learning_rate'] = 200
        d['manifold']['args']['n_iter'] = 5000
        d['manifold']['args']['n_iter_without_progress'] = 300,
        d['manifold']['args']['min_grad_norm'] = 1e-7,
        # d['manifold']['args']['metric']='euclidean',
        d['manifold']['args']['init'] = 'pca',
        d['manifold']['args']['verbose'] = 0,
        d['manifold']['args']['random_state'] = None
        # method='barnes_hut'
        # angle=0.5

        d['neighbors'] = OrderedDict()
        d['neighbors']['type'] = 'ball_tree'
        d['neighbors']['kNN'] = 4
        d['neighbors']['args'] = OrderedDict()
        d['neighbors']['args']['leaf_size'] = 40
        d['neighbors']['args']['metric'] = 'minkowski'

        d['cluster'] = OrderedDict()
        d['cluster']['type'] = 'dbscan'
        d['cluster']['args'] = OrderedDict()
        d['cluster']['args']['eps'] = OrderedDict()
        d['cluster']['args']['eps']['NN'] = 3
        d['cluster']['args']['eps']['percentile'] = .99
        d['cluster']['args']['min_samples'] = 10
        d['cluster']['args']['metric'] = 'euclidean'
        d['cluster']['args']['metric_params'] = None
        d['cluster']['args']['algorithm'] = 'auto'
        d['cluster']['args']['leaf_size'] = 30
        d['cluster']['args']['p'] = None
        # add conditional here to check for custom clustering params in the config
        # pseudo: if 'clustering_params' in self.configuration:
        # pseudo:    for k,v in self.configuration['clustering_params'].items():
        # pseudo:        d[k] = v
        # not yet implemented in configuration
        return d


if __name__ == "__main__":
    o = PyposmatClusterSampler()
    # read the configuration file
    o.read_configuration_file(
            filename=pyposmat_filename_in
    )
    # read the pyposmat datafile
    o.configure_pyposmat_datafile_in(
        filename=pyposmat_datafile_in
    )
    i_iteration = 0
    o.run_simulations(i_iteration=0)

