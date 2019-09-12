"""
pypospack.dft

This is module for some functional utilies for DFT simulations
"""

import copy
import numpy as np
from collections import OrderedDict
import pypospack.crystal as crystal

def get_kpoint_mesh(simulation_cell,linear_kpoint_density):
    """

    This function determines what the kpoint mesh should be given a linear
    kpoint density.  This is inspired from the kpoint generation methodology
    of VASP.

    Args:
        simulation_cell(pyflamestk.crystal.SimulationCell):
        linear_kpoint_density(str)

    Refs:
        https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
    """

    # assert isinstance(simulation_cell,crystal.SimulationCell)
    b1_length = np.dot(simulation_cell.b1,simulation_cell.b1)**0.5
    b2_length = np.dot(simulation_cell.b2,simulation_cell.b2)**0.5
    b3_length = np.dot(simulation_cell.b3,simulation_cell.b3)**0.5

    #linear_kpoint_density = kpoints_1/b1_length -->i
    kpoints_1 = linear_kpoint_density*b1_length
    kpoints_2 = linear_kpoint_density*b2_length
    kpoints_3 = linear_kpoint_density*b3_length

    kpoints_1_round = int(kpoints_1+0.5)
    kpoints_2_round = int(kpoints_2+0.5)
    kpoints_3_round = int(kpoints_3+0.5)

    return [kpoints_1_round,kpoints_2_round,kpoints_3_round]

def determine_kpoint_meshes(simulation_cell,
                            rho_min=1,
                            rho_max=10,
                            d_rho=0.1,
                            kpoint_min=3,
                            kpoint_max=15):
    """

   Args:
        simulation_cell(pypospack.crystal.SimulationCell):  The simulation
            cell in which to determine the kpoint density.  Any class which
            subclasses pypospack.crystal.SimulationCell is usable.
        rho_min(float): lowest linear kpoint density to start. Default is 1.
        rho_max(float): largest linear kpoint density to end.  Default is 10.
        d_rho(float): the size of the step of the linear kpoint density. Default
            is 0.1
        kpoint_min(int): Minimum number of kpoints in a linear direction.
            Default is 3, which is appropriate for approximately cubic
            crystals.
        kpoint_max(int): Maximum number of kpoints is a linear direction.
            Default is 15, which is appropriate for approximately cubic
            crystals.

    Returns:
        OrderedDict:
            key(str): 'kp_{k1}_{k2}_{k3}'
            values(list of int): [{k1},{k2},{k3}]
    """

    def kpoints2key(kpm):
        return "kp_{}_{}_{}".format(kpm[0], kpm[1], kpm[2])

    assert isinstance(simulation_cell,crystal.SimulationCell)
    # initialize some variables
    old_kp_mesh = None
    new_kp_mesh = None
    kpoint_meshes = OrderedDict()

    # search for appropriate kpoint meshes
    for rho in np.arange(rho_min,rho_max,d_rho):
        old_kp_mesh = copy.copy(new_kp_mesh)
        new_kp_mesh = get_kpoint_mesh(simulation_cell,rho)

        if old_kp_mesh is None:
            key = kpoints2key(new_kp_mesh)
            kpoint_meshes[key] = list(new_kp_mesh)
        else:
            is_kpoints_unchanged = all([
                old_kp_mesh[0] == new_kp_mesh[0],
                old_kp_mesh[1] == new_kp_mesh[1],
                old_kp_mesh[2] == new_kp_mesh[2]])
            if not is_kpoints_unchanged:
                key = kpoints2key(new_kp_mesh)
                kpoint_meshes[key] = list(new_kp_mesh)
    return copy.deepcopy(kpoint_meshes)
