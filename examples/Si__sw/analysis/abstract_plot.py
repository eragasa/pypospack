import matplotlib.pyplot as plt
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

class PyposmatAbstractPlot(object):

    def __init__(self,config=None,data=None):

        self.configuration = None
        self.data = None
        self.fig = None
        self.ax = None

        self.initialize_configuration(config=config)
        self.initialize_data(data=data)

    def create_subplots(self,nrows=1,ncols=1,
            sharex=False,sharey=False,
            squeeze=True,
            subplot_kw=None,gridspec_kw=None):
        """ create subplots 
        
        has the same arguments as the matplotlib.pyplot.subplots command
        """
        plt.close()
        self.fig, self.ax = plt.subplots(
                nrows = nrows,
                ncols = ncols,
                sharex = sharex,
                sharey = sharey,
                squeeze = True,
                subplot_kw = subplot_kw,
                gridspec_kw = gridspec_kw)

    def initialize_configuration(self,config):
        assert isinstance(config,str) \
                or isinstance(config,PyposmatConfigurationFile) \
                or config is None

        if isinstance(config,str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=config)
        elif isinstance(config,PyposmatConfigurationFile):
            self.configuration = config
        elif config is None:
            self.configuration = None
        else:
            m = 'config arguement must either be a path string of a PyposmatConfigurationFile object'
            raise TypeError(m)
        
    def initialize_data(self,data):
        assert isinstance(data,str) \
                or isinstance(data,PyposmatDataFile) \
                or data is None

        if isinstance(data,str):
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        elif isinstance(data,PyposmatDataFile):
            self.data = data
        elif data is None:
            self.data = None
        else:
            m = 'data argument must either be path string or a PyposmatDataFile object'
            raise TypeError(m)


    def savefig(self,filename,dpi=1200,transparent=True):
        
        self.fig.savefig(
                fname=filename,
                dpi=dpi,
                transparent=True)



