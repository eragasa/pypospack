import math
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pypospack.pareto as pareto
n_simulations = 0


def kurasawe(x):
    global n_simulations
    n_simulations+=1
    return kurasawe_f1(x), kurasawe_f2(x)

def kurasawe_f1(x):
    n = len(x)

    f1 = 0
    for i in range(n-1):
        f1 += -10*math.exp(-0.2*(math.sqrt(x[i]**2 + x[i+1]**2)))
    return f1

def kurasawe_f2(x):
    n = len(x)

    f2 = 0
    for i in range(n):
        f2 += abs(x[i])**0.8+5*math.sin(x[i]**3)

    return f2

def kurasawe_cost_function(x,w1,w2):
    global n_simulations
    weights = (w1,w2)

    obj_functions = [kurasawe_f1,kurasawe_f2]
    C = 0
    for i in range(len(obj_functions)):
        C += weights[i] * obj_functions[i](x)

    n_simulations = n_simulations + 1
    return C


class ParetoProblem(object):
    pass

class KurasaweParetoProblem(ParetoProblem):
    def __init__(self):
        self.design_variables = (None,None,None)

        self.design_variable_ranges = (
            (-5,10) for k in len(self.design_variables)
        )
        self.objective_functions = (kurasawe_f1,kurasawe_f2)

        self.cost_function = kurasawe_cost_function

class ParetoSurfaceEstimatorConjugateGradientMethod():
    def __init__(self,pareto_problem):
        assert isinstance(pareto_problem,ParetoProblem)

def costfunction(x,w1,w2):
    global n_simulations
    weights = (w1,w2)

    obj_functions = [kurasawe_f1,kurasawe_f2]
    C = 0
    for i in range(len(obj_functions)):
        C += weights[i] * obj_functions[i](x)

    n_simulations = n_simulations + 1
    return C

def pareto_surface_from_conjugate_gradient_method():
    from scipy.optimize import minimize

    weights_max = 1.00
    weights_N = 10

    n_design_variables = 3
    n_obj_functions = 2

    weights = []
    weights_range = np.array(
            weights_max/weights_N*np.linspace(1,weights_N,weights_N)
        ).tolist()

    column_names = ['id'] \
            + ['x{}'.format(i+1) for i in range(n_design_variables)] \
            + ['f{}'.format(i+1) for i in range(n_obj_functions)]
    data = []
    for w1 in weights_range:
        for w2 in weights_range:
            w_norm = np.sqrt(w1**2 + w2**2)
            weights.append((w1,w2))

    x1_min = x2_min = x3_min = -5
    x1_max = x2_max = x3_max = 10
    x1_N = x2_N = x3_N = 5

    i=0
    for x1 in np.linspace(x1_min,x1_max,x1_N).tolist():
        for x2 in np.linspace(x2_min,x2_max,x2_N).tolist():
            for x3 in np.linspace(x3_min,x2_max,x3_N).tolist():
                print("\tx=({x1},{x2},{x3})".format(x1=x1,x2=x2,x3=x3))
                for w in weights:
                    optimize_result = minimize(fun=costfunction,
                             x0 = np.array([x1,x2,x3]),
                             args=(w),
                             method='L-BFGS-B'
                    )
                    xopt = optimize_result.x
                    f1 = kurasawe_f1(xopt)
                    f2 = kurasawe_f2(xopt)
                    results = [i]+xopt.tolist()+[f1,f2]
                    data.append(results)
                    i += 1
    df = pd.DataFrame(data, columns=column_names)
    return df

def pareto_surface_from_monte_carlo_sampling():
    from scipy.stats import uniform

    n_design_variables = 3
    n_obj_functions = 2
    N_simulations = 100000

    x_sampler = [uniform(-5,10) for i in range(n_design_variables)]

    column_names = ['id'] \
            + ['x{}'.format(i+1) for i in range(n_design_variables)] \
            + ['f{}'.format(i+1) for i in range(n_obj_functions)]
    data = []
    for i in range(N_simulations):
        x = [u.rvs() for u in x_sampler]
        f1,f2 = kurasawe(x)
        row = [i] + x + [f1,f2]
        data.append([i] + x + [f1,f2])

    # return the results as a pandas.DataFrame
    df = pd.DataFrame(data, columns=column_names)
    return df

def pareto_surface_from_kde_sampling():
    from scipy.stats import uniform, gaussian_kde

    n_design_variables = 3
    n_obj_functions = 2
    N_iterations = 10
    N_samples_per_iteration = 10000

    column_names = ['id'] \
            + ['x{}'.format(i+1) for i in range(n_design_variables)] \
            + ['f{}'.format(i+1) for i in range(n_obj_functions)]

    data = []
    df = None
    for i_iteration in range(N_iterations):
        if i_iteration == 0:
            x_sampler = [uniform(-5,10) for i in range(n_design_variables)]
            for i in range(N_samples_per_iteration):
                x = [u.rvs() for u in x_sampler]
                f1,f2 = kurasawe(x)
                row = [i] + x + [f1,f2]
                data.append([i] + x + [f1,f2])
                df = pd.DataFrame(data, columns=column_names)

        else:
            pareto_rows = pareto.pareto(df[['f1','f2']].values.tolist())
            print("i={},n_pareto={}".format(i_iteration,len(pareto_rows)))
            print(df.loc[pareto_rows,['x1','x2','x3']].T)
            kde = gaussian_kde(df.loc[pareto_rows,['x1','x2','x3']].T)
            data = []
            for i in range(N_samples_per_iteration):
                x = kde.resample(1)

                f1,f2 = kurasawe([x[i] for i in range(3)])
                row = [i] + x + [f1,f2]
                data.append([i] + [x[i][0] for i in range(3)] + [f1,f2])
            df = pd.DataFrame(data, columns=column_names)
    return df

fig, ax = plt.subplots(1,3,figsize=(5,5/3))
if True:
    n_simulations = 0
    cg_df = pareto_surface_from_conjugate_gradient_method()
    cg_pareto_rows = pareto.pareto(cg_df[['f1','f2']].values.tolist())
    cg_pareto_df = cg_df.loc[cg_pareto_rows]
    cg_pareto_df = cg_pareto_df.sort_values('f1')
    ax[0].scatter(
        cg_df['f1'],
        cg_df['f2'],
        c='black',marker=".",s=1
    )
    ax[0].plot(
        cg_pareto_df['f1'],
        cg_pareto_df['f2'],
        c='red',linewidth=2,marker=".",ms=2
    )
    print(n_simulations)
    print(len(cg_pareto_rows))

if True:
    n_simulations = 0
    mc_df = pareto_surface_from_monte_carlo_sampling()
    #mc_df.sort_values(by=['f1'], ascending=True,inplace=True, kind='merges')
    mc_pareto_rows = pareto.pareto(mc_df[['f1','f2']].values.tolist())
    mc_pareto_df = mc_df.loc[mc_pareto_rows]
    mc_pareto_df = mc_pareto_df.sort_values('f1')
    ax[1].scatter(
        mc_df['f1'],
        mc_df['f2'],
        c='black',marker=".",s=1
    )
    ax[1].plot(
        mc_pareto_df['f1'],
        mc_pareto_df['f2'],
        #mc_df.loc[mc_pareto_rows,'f1'],
        #mc_df.loc[mc_pareto_rows,'f2'],
        c='red',linewidth=2,marker=".",ms=2
    )


    print(n_simulations)
    print(len(mc_pareto_rows))
if True:
    n_simulations = 0
    kde_df = pareto_surface_from_kde_sampling()
    kde_pareto_rows = pareto.pareto(kde_df[['f1','f2']].values.tolist())
    kde_pareto_df = kde_df.loc[kde_pareto_rows]
    kde_pareto_df = kde_pareto_df.sort_values('f1')
    ax[2].scatter(
        kde_df['f1'],
        kde_df['f2'],
        c='black',marker=".",s=1
    )
    ax[2].plot(
        kde_pareto_df['f1'],
        kde_pareto_df['f2'],
        #mc_df.loc[mc_pareto_rows,'f1'],
        #mc_df.loc[mc_pareto_rows,'f2'],
        c='red',linewidth=2,marker=".",ms=2
    )
    print(n_simulations)
    print(len(kde_pareto_rows))

for i in range(3):
    ax[i].set_xlim(-20,0)
    ax[i].set_ylim(-12.5,7.5)
    ax[i].set_xlabel("f_1(x)")
    ax[i].set_ylabel("f_2(x)")
    ax[i].set(adjustable='box-forced',aspect='equal')
plt.tight_layout()
fig.savefig("comparison.eps")
plt.show()
