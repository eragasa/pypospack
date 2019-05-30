import numpy as np
import scipy.optimize

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_cg.html

class ObjectiveFunction(object):
    def __init__(self):
        pass

    def evaluate(self,x,*args):
        u, v = x
        a, b, c, d, e, f = args
        return a*u**2 + b*u*v + c*v**2 + d*u + e*v + f

    def callback(self,xk):
        print(xk)

class MultiObjectiveFunction(object):
    " multi-objective optimization problem "
    def __init__(self):
        self.objective_functions = []

    def evaluate(self,theta):
        return [f(theta) for k,f in self.objective_functions.items]
    
    def callback(self,xk):
        pass

class BinhAndKornFunction(MultiObjectiveFunction):
    def __init__(self):
        MultiObjectiveFunction.__init__(self)
        self.objective_functions = OrderedDict
        self.objective_functions['f1'].append(self.func_f1)
        self.objective_functions['f2'].append(self.func_f1)

    def func_f1(self,theta):
        x = theta[0]
        y = theta[1]
        return 4*x**2 + 4*y**2

    def func_f2(self,theta):
        x = theta[0]
        y = theta[1]
        return (x-5)**2 + (y-5)**2
        
class SingleObjectiveOptimization(object):
    pass

class MultiObjectiveOptimization(object):
    pass

class WeightedCostFunction(ObjectiveFunction):
    def __init__(self,moo_problem,weights):
        assert isinstance(moo_problem,MultiObjectiveFunction)
        self.objective_function = moo_problem

        assert len(weights) == len(moo_problem.objective_functions)
        self.weights = weights

    def evaluate(self,theta):
        obj_func_vals = self.objective_function.evaluate(theta)
        
        weighted_cost = 0
        try:
            for i,w in enumerate(self.weights):
                weighted_cost += w * obj_func_vals[i]
        except:
            print("weights:{}".format(self.weights))
            raise

        values_to_str = [str(v) for v in theta.tolist()+obj_func_vals]
        print(",".join(values_to_str))
        return weighted_cost

class ParetoOptimization(MultiObjectiveOptimization):
    def __init__(self,moo_problem):
        assert isinstance(moo_problem,MultiObjectiveFunction)
        self.objective_function = moo_problem

    def minimize(self,theta0):
        pass

class ConjugateGradientOptimization(SingleObjectiveOptimization):
    def __init__(self,obj_function):
        SingleObjectiveOptimization.__init__(self)
        self.objective_function = obj_function

    def minimize(self,theta0):
        result = scipy.optimize.fmin_cg(\
                self.objective_function.evaluate,
                theta0)

if __name__ == '__main__':
    w = [0.5, 0.5]
    theta0 = [1,1]
    bkf = BinhAndKornFunction()
    w_bkf = WeightedCostFunction(bkf,w)
    cg_w_bkf = ConjugateGradientOptimization(w_bkf)
    cg_w_bkf.minimize(theta0)
