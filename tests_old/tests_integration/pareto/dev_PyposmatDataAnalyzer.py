# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, copy
import numpy as np
from collections import OrderedDict()
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatDataAnalyzer
import pypospack.pareto as pareto



class ParetoProblem(object):
    pass
class ParetoConesProblem(ParetoProblem):

    def __init__(self):
        ParetoProblem.__init__(self)
        self.objective_function = None
        self.variables = None

        self.add-
    def add_variables(self)
        self.variables = OrderedDict()
        self.variables['r'] = None
        self.variables['h'] = None

    def add_variable_constraints(self):
        self.variable_constraints = OrderedDict()
        self.variable_constraints['r<10'] = 'r<10'
        self.variable_constraints['r>0'] = 'r>10'
        self.variable_constraints['h>0'] = 'h>0'
        self.variable_constraints['h<20'] = 'h<20'

    def add_objective_functions(self):
        self.objective_functions=OrderedDict()
        self.objective_function['lateral_surface_area'] \
            = lateral_surface_area
        self.objective_function['total_surface_area'] \
            = total_surface_area

    def add_constraints(self):
        self.constraints=OrderedDict()
        self.constraints['volume']

    def base_surface_area(r,s):
        return np.pi*r**2
    def lateral_surface_area(r,s):
        return np.pi*r*s
    def total_surface_area(r,s):
        B = base_surface_area(r,s)
        S = lateral_surface_area(r,s)
        T = total_surface_area(r,s)
        return
    def volume(r,s):
        V = np.pi/3*(r**2)*h
        return V
    def slant_height(r,h):
        return np.sqrt(r*r+h*h)

if __name__ == "__main__":

    data_directory = 'data'
    pyposmat_data_filename = 'pypospack.results.out'
    pyposmat_configuration_filename = 'pypospack.config.in'
    data_analyzer = PyposmatDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=pyposmat_configuration_filename)
    data_analyzer.read_data_file(
            filename=os.path.join(
                data_directory,
                pyposmat_data_filename))
    data_analyzer.calculate_pareto_set()
    data_analyzer.write_kde_file(
            filename=os.path.join(
                data_directory,
                'pyposmat.kde.out'))

    exit()
