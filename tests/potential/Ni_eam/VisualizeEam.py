import matplotlib.pyplot as plt
import pypospack.potential as potential
import numpy as np

class EamPotentialFigure(object):
    def __init__(self):
        self.figure = plt.figure()
        self.plot_pair = plt.subplot(3,1,1)
        self.plot_density = plt.subplot(3,1,2)
        self.plot_embedding = plt.subplot(3,1,3)

        self.pair_func = None
        self.embed_func = None
        self.elec_dens_func = None

    def plot_eam_potential(self,o_eam_potential, parameters):

    def create_figure(self):

class EamPairPotentialPlot(object):
    pass

class EamDensityFunctionPlot(object):
    pass

class EamEmbeddingFunctionPlot(object):
    pass
