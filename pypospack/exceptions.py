
"""
pypospack exception classes
"""

class BaseException(Exception):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.msg = args[0]

    def __str__(self):
        return self.msg

    def explain(self):
        return "{}: {}".format(
            self.__class__.__name__,
            self.msg)

class BadPreprocessorTypeException(BaseException): pass

class BadParameterException(BaseException): pass

class BadManifoldTypeException(BaseException): pass

class BadClusterTypeException(BaseException): pass

class BadNearestNeighborTypeException(BaseException): pass

class LammpsSimulationError(BaseException): pass

class PyposmatError(BaseException): pass
