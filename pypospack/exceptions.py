
"""
pypospack exception classes
"""

class BaseException(Exception):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        try:
            self.msg = args[0]
        except IndexError as e:
            self.msg = ""

    def __str__(self):
        return self.msg

    def explain(self):
        return "{}: {}".format(
            self.__class__.__name__,
            self.msg)

#TODO: rneame *Exception names to *Error names

class BadPreprocessorTypeException(BaseException): pass

class PyposmatBadParameterError(BaseException): pass
class BadParameterException(PyposmatBadParameterError): 
    pass
# old class saved here just in case
#class BadParameterException(Exception):
#    def __init__(self,
#            code,
#            parameter_name,
#            parameter_value,
#            parameters):
#        self.code = code
#        self.parameter_name = parameter_name
#        self.parameter_value = parameter_value
#        self.parameters = parameters

class BadManifoldTypeException(BaseException): pass

class BadClusterTypeException(BaseException): pass

class BadNearestNeighborTypeException(BaseException): pass

class LammpsSimulationError(BaseException): pass

class PyposmatError(BaseException): pass

class PypospackTaskManagerError(BaseException): pass

class PypospackBadKdeBandwidthType(BaseException): pass

class PyposmatUnknownDistributionType(BaseException): pass
