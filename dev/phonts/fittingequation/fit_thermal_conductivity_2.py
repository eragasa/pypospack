from collections import OrderedDict

import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import metrics

filename = "ThermalConductivity_temperature.dat"

def read_thermal_conductivity_data(path):
    """imports the thermal conductivity data"""
    with open(path,'r') as f:
        lines = [k.strip() for k in f.readlines()]

    rows = []
    for i, line in enumerate(lines):
        line_args = [k.strip() for k in line.split()]
        if i == 0:
            column_names = line_args
        else:
            row = [float(k) for k in line_args]
            rows.append(row)
    df = pd.DataFrame(data=rows,
                      columns=column_names)

    return df

def thermal_conductivity_formula(x, k0, alpha, beta):
    temperature = x[:,0]
    pressure = x[:,1]
    return k0 * (1+beta*pressure) / (1+alpha*temperature)



if __name__ == "__main__":
    df = read_thermal_conductivity_data(path=filename)

    df['1/T'] = 1/df['T']
    pressures = list(set(df['P'].tolist()))
    temperatures = list(set(df['T'].tolist()))

    pressures.sort()
    temperatures.sort()

    kxys = []
    for i in range(1,4):
        kxy = 'k{}{}'.format(i,i)
        df['1/{}'.format(kxy)] = 1/df[kxy]
        kxys.append(kxy)

    alphas = OrderedDict()
    betas = OrderedDict()

    beta_regression_columns = ['kxy','P', 'b1', 'b2' 'beta']
    beta_regression_data = []

    fig, axes = plt.subplots(1,3)
    for i, kxy in enumerate(kxys):
        betas[kxy] = OrderedDict()
        print(kxy)
        for P in pressures:
            regressor = LinearRegression()

            df_P = df.loc[df['P'] == P]
            X = df_P['T'].values.reshape(-1,1)
            Y = 1/df_P[kxy].values.reshape(-1,1)
            regressor.fit(X,Y)

            b1 = float(regressor.intercept_)
            b2 = float(regressor.coef_)
            beta = b2/b1

            beta_regression_data.append([kxy, P, b1, b2, beta])
            
            betas[kxy][P] = beta

            axes[i].scatter(X,Y)
            axes[i].plot(np.array(temperatures), 
                                  b1 + b2 * np.array(temperatures))
        axes[i].set_xlim(300,900)
        axes[i].set_xlabel('1/T [K]')
        axes[i].set_ylabel(kxy)
        axes[i].axis('auto')

    plt.tight_layout()
    plt.show()
    
    beta_regression_df = pd.DataFrame(data=beta_regression_data,
                                      columns=beta_regression_columns)


    alpha_regression_columns = ['kxy', 'T', 'a1', 'a2', 'alpha']
    alpha_regression_data = []

    fig, axes = plt.subplots(1,3)
    for i,kxy in enumerate(kxys):
        alphas[kxy] = OrderedDict()
        for T in temperatures:
            regressor = LinearRegression()

            df_T = df.loc[df['T'] == T]
            X = df_T['P'].values.reshape(-1,1)
            Y = df_T[kxy].values.reshape(-1,1)
            regressor.fit(X,Y)

            a1 = float(regressor.intercept_)
            a2 = float(regressor.coef_)
            alpha = a2/a1
            alphas[kxy][T] = alpha

            alpha_regression_data.append([kxy, T, a1, a2, alpha])

            axes[i].scatter(X,Y)
        axes[i].set_xlim(0,4)
        axes[i].set_xlabel('P [GPa]')
        axes[i].set_ylabel(kxy)
        axes[i].axis('auto')
    plt.tight_layout()
    plt.draw()
    plt.show()

    alpha_regression_df = pd.DataFrame(data=alpha_regression_data,
                                       columns=alpha_regression_columns)

    k0 = np.zeros(shape=(len(temperatures),len(pressures)))
    coeff = np.zeros(shape=(len(temperatures), len(pressures)))

    def surface_plot (matrix, **kwargs):
        # acquire the cartesian coordinate matrices from the matrix
        # x is cols, y is rows
        (x, y) = np.meshgrid(np.arange(matrix.shape[0]), np.arange(matrix.shape[1]))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(x, y, matrix, **kwargs)
        return (fig, ax, surf)

    def get_k0(T,P,kxy):
        T = temperatures[m]
        P = pressures[n]
        alpha = alphas[kxy][T]
        beta = betas[kxy][P]
        k = df[(df['P'] == P) & (df['T'] == T)]
        k = float(k[kxy])
        k0[iT, iP] = (1+alpha*T)/(1+beta*P)*k
        coeff[iT, iP] = (1+beta*P)/(1+alpha*T)
        return(k0)

    for i,kxy in enumerate(kxys):
        for iT, T in enumerate(temperatures):
            for iP, P in enumerate(pressures):
                alpha = alphas[kxy][T]
                beta = betas[kxy][P]
                k = df[(df['P'] == P) & (df['T'] == T)]
                k = float(k[kxy])
                k0[iT, iP] = (1+alpha*T)/(1+beta*P)*k
                coeff[iT, iP] = (1+beta*P)/(1+alpha*T)
        T, P = np.meshgrid(temperatures, pressures)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(T,P,k0.T)
        ax.set_xlabel('T [K]')
        ax.set_ylabel('P [GPa]')
        ax.set_zlabel('k0')
        plt.show()
        print(k0)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(T,P,coeff.T)
        ax.set_xlabel('T [K]')
        ax.set_ylabel('P [GPa]')
        ax.set_zlabel('coeff')

        plt.show()
        print(coeff)
