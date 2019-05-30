import pytest
import optimize

if __name__ == '__main__':
    weights = [0.5,0.5]
    theta0 = [1,1]
    func_binhkorn = optimize.BinhAndKornFunction()
    func_wcf_binhkorn = optimize.WeightedCostFunction(
            func_binhkorn, weights)
    cg_opt = optimize.ConjugateGradientOptimization(
            obj_function = func_wcf_binhkorn)
    cg_opt.minimize(theta0)
