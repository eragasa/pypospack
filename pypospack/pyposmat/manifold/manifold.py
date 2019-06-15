# -*- coding: utf-8 -*-
"""Implementation of abstract manifold class

This module provides the abstract class implementation for manifold learning
this wrapper is primarily a wrapper from the scikit-learn module.
"""
from sklearn import preprocessing

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class Manifold():
    """abstract manifold class

    Args:
        pyposmat_configuration (:obj:`str`,:obj:`PyposmatConfigurationFile`):
            Either an instance of a pypospack configuration file or the path to
            an appropriate file
        pyposmat_data (:obj:`str`,:obj:`PyposmatDataFile`): Either an instance
            of a pyposmat data file or the path to an appropriate file
        manifold_configuration (:obj:`OrderedDict`,optional): appropriate
            configurations for the manifold. If set to None, the implemented
            Manifold object should provide default settings.

    Attributes:
        configuration (PyposmatConfigurationFile): an instance of the
            configuration object.
        data (PyposmatDataFile): an instance of the data file object
        manifold_configuration (OrderedDict): contains the arguments for the
            manifold learning technique
    """

    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 manifold_configuration=None):

        # define attributes
        self.configuration = None
        self.data = None
        self.manifold_configuration = None

        # initialization segment
        self.initialize_configuration(configuration=pyposmat_configuration)
        self.initialize_data(data=pyposmat_data)
        self.initialize_manifold_configuration(configuration=manifold_configuration)

    def initialize_configuration(self, configuration):
        """ initialize pyposmat configuration

        For the convenience purposes, definitions of QOI names and parameter
        names are taken from the configuration file object.

        Args:
            configuration (:obj:`str`,:obj:`PyposmatConfigurationFile`): Either
                an instance of the :obj:`PyposmatConfigurationfile` or a path
                to a valid file
        """

        if isinstance(configuration, PyposmatConfigurationFile):
            self.configuration = configuration
        elif isinstance(configuration, str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        else:
            msg = ('configuration must be either a path or '
                   'PyposmatConfigurationFile')
            raise TypeError(msg)

    def initialize_data(self, data):
        """ initialize pyposmat data

        For convenience purposes, the datafile from the appropriate iteration
        from a pypospack simulation.

        Args:
            configuration (:obj:`str`,:obj:`PyposmatDataFile`): Either an
                an instance of the :obj:`PyposmatDataFile` or a path to a valid
                file.
        """

        if isinstance(data, PyposmatDataFile):
            self.data = data
        elif isinstance(data, str):
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        else:
            msg = 'data must either be a path or a PyposmatDataFile'
            raise TypeError(msg)

    def initialize_manifold_configuration(self, configuration=None):
        """ initialize manifold configuration

        The manifold configuraiton is stored internally as an :obj:`OrderedDict`
        as the `manifold_configuration` attribute.  Classes which inherit from
        this class should override this implementation.

        Args:
            configuration(None,OrderedDict): a dictionary type object which
                contains all the information of the manifold configruation
        """

        if configuration is None:
            self.manifold_configuration = None
        else:
            raise NotImplementedError

    def learn_manifold(self, names, scaling_type):
        """ learn the manifold

        Args:
            names (list of str, 'qoi','free_parameters','all'): list for the
                column names, which should be included in the manifold learning
                analysis.  'qoi' uses all the qoi names as defined in the
                configuration attributes. 'free_parameter' uses the free
                parameter names determine from the configuration attribute.
                'all' combines the qoi names and the free parameter names.
                Alternatively, a list of column names can be passed
            scaling_type (str): this is dependent upon manifold learning type.
                Different manifold learning types perform better due to
                guidelines.
        """

        raise NotImplementedError()

    def scale_data(self, X, scaling_type): # pylint: disable=invalid-name
        """ scales data

        Args:
            X (numpy.ndarray,pandas.DataFrame): the data
            scaling_type (str): the scaling type
        """

        if scaling_type == 'standard':
            X_scaled = preprocessing.scale(X) # pylint: disable=invalid-name
        elif scaling_type == 'none':
            X_scaled = X # pylint: disable=invalid-name
        else:
            raise ValueError('unknown scaling type')

        return X_scaled
if __name__ == "__main__":
    import os
    import pypospack.utils
    # pylint: disable=invalid-name

    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

    config_fn = os.path.join(
        pypospack_root_dir,
        'data', 'Si__sw__data', 'pareto_optimization_unconstrained',
        'pyposmat.config.in')
    data_fn = os.path.join(
        pypospack_root_dir,
        'data', 'Si__sw__data', 'pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    o = Manifold(pyposmat_configuration=o_config,pyposmat_data=o_data)
