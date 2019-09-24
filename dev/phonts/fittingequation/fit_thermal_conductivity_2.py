import numpy as np
import pandas as pd
from scipy import optimize
from sklearn.linear_model import LinearRegression
filename = "ThermalConductivity_temperature.dat"

def read_thermal_conductivity_data(path):
    with open(path,'r') as f:
        lines = [k.strip() for k in f.readlines()]

    rows = []

    for i line in enumerate(lines):
        line_args = [k.strip() for k in line.split()]
        if i == 0:
            column_names = line_args
        else:
            row = [float(k) for k in line_args]

    df = pd.DataFrame(data=rows,
                      columns=column_names)
    
    return df

def thermal_conductivity_formula(x, k0, alpha, beta):
    temperature = x[:,0]
    pressure = x[:,1]
    return k0 * (1+beta*pressure) / (1+alpha*temperature)

def log_thermal_conductivity_formala(x, k0, alpha, beta):
    ln_temperature = np.log(x[:,0])
    ln_pressure = np.log(x[:,1])

    return np.log(k0) + np.log(1+beta*pressure) - np.log(1+alpha*temperature)


for k in ['k11','k22','k33']:
    print(80*'-')
    print(k)
    y = df[k].values.tolist()
    x = df[['T','P']].values.tolist()

    popt, pcov = optimize.curve_fit(thermal_conductivity_formula,x,y)

    print(popt)
    print(pcov)
