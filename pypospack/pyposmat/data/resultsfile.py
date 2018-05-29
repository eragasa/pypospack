import copy
import yaml
from collections import OrderedDict
from pypospack.io.filesystem import OrderedDictYAMLLoader

class PyposmatResultsFile(object):

    def __init__(self):
        self.filename_in = None
        self.filename_out = None
        self.results = None

        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None

    
