from pypospack.pareto import pareto
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pareto_visualization import make_pareto_plot_2d
if __name__ == "__main__":
    x = np.random.normal(0,1,1000)
    y = np.random.normal(0,1,1000)

    df = pd.DataFrame()
    df['x'] = x
    df['y'] = y

    pareto_idx = pareto(df.values.tolist())

    make_pareto_plot_2d(df,pareto_idx)
