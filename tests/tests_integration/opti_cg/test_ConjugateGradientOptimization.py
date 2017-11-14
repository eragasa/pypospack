import pytest
import optimize

w = [0.5,0.5]
theta0 = [1,1]

def test_BinhAndKornFunction__init():
    func_binhkorn = optimize.BinhAndKornFunction()

def test_WeightCostFunction__init():
    weights = [0.5,0.5]
    func_binhkorn = optimize.BinhAndKornFunction()
    wcf_bkf = optimize.WeightedCostFunction(
            func_binhkorn,
            weights)

def test_ConjugateGradientOptimization__init():
    weights = [0.5,0.5]
    theta0 = [1,1]
    func_binhkorn = optimize.BinhAndKornFunction()
    func_wcf_binhkorn = optimize.WeightedCostFunction(
            func_binhkorn, weights)
    cg_opt = optimize.ConjugateGradientOptimization(
            obj_function = func_wcf_binhkorn)
    cg_opt.minimize=theta0
