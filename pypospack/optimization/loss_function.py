
class LossFunction(object):
    """ loss function base object

    a loss function or cost function is a function that maps an event or values of one or more variables onto a real number intuitively representing some "cost" associated with the event.
    """
    def evaluate(self):
        assert isinstance(qoi_predicted,float)
        assert isinstance(qoi_target,float)

class AbsoluteLossFunction(LossFunction):

    def evaluate(self,qoi_predicted,qoi_target):
        assert isinstance(qoi_predicted,float)
        assert isinstance(qoi_target,float)

        return abs(qoi_predicted-qoi_target)

class QuadraticLossFunction(LossFunction):

    def evaluate(self,qoi_predicted,qoi_target):
        assert isinstance(qoi_predicted,float)
        assert isinstance(qoi_target,float)
        
        return (qoi_predicted-qoi_target)**2.

class WeightCostFunction(LossFunction):

    def add_loss_function(self,qoi_name,weight,loss_function):
        assert isinstance(weight,float)
        assert isinstance(loss_function,LossFunction)

        self.loss_functions[qoi_name] = OrderedDict()
        self.loss_functions[qoi_name]['weight'] = weight
        self.loss_functions[qoi_name]['function'] = loss_function

    def evaluate_qois(self,parameters):

    def evaluate(self):

        cost_function = 0.
        for loss_function in self.loss_functions:
            cost_function += loss_function['weight'] * loss_function['function'].evaluate()

        return cost_function


