import os,copy
from collection import OrderedDict
from pypospack.qoi import QoiManager
from pypospack.task import TaskManager
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PYposmatDataFile

class SlurmSubmissionEngine(object):


