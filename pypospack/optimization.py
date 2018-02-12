""" classes and functions for optimization techniques """

class MultiObjectiveProblem(object):
    """ A multiple objective function class

    A multiple objective problem is defined as a problem which requires the
    simultaneous minimization of multiple objectives at the same time.

    Attributes:
        obj_func (list): a list of objective functions
    """
    def __init__(self):
        self.obj_func

class ObjectiveFunction(object):


class MultiobjectiveParameterOptimization(object):
    pass

class ParetoParameterOptimzation(MultiobjectiveParameterOptimization):
    def __init__(self):
        pass

class ConjugateGradientParameterOptimization(MultiobjectiveParameterOptimization):
    pass
