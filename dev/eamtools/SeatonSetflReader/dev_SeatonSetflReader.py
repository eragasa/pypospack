#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = "Eugene Ragasa"
__copyright__ = "Copyright 2019, Eugene J. Ragasa"
__licence__ = "MIT"
__contributors__ = ["Eugene J Ragasa, R. Seaton Ullberg, Simon Phillpot"]
__maintainer__ = "Eugene J. Ragasa"

# testing library
import pytest

# standard libraries
import os
from collections import OrderedDict

# third part libraries
import numpy as np

# local source
import pypospack.utils
from pypospack.eamtools import (SeatonSetflReader,
                                get_setfl_pair_order)

if False:
    def create_r(r_max,n):
        r = r_max/n * np.linspace(1,n,n)
        return r

    def get_setfl_pair_order(symbols):
        pairs = []
        for i,s_i in enumerate(symbols):
            for j,s_j in enumerate(symbols):
                if i >= j:
                    pairs.append([s_i,s_j])
        return pairs

def plot_eam_potential_from_setflfile(filename):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3,1)

    setfl = SeatonSetflReader(path=setfl_filename)
    setfl.read()

    # plot the pair potential function
    for pair in setfl.element_pairs:
        pair_x = create_r(setfl.n_r*setfl.d_r,setfl.n_r)
        pair_y = setfl.pair_function(pair)/pair_x

        ax[0].plot(pair_x,pair_y,label=pair)
    ax[0].legend()
    
    # plot the density functions
    for s in setfl.elements:
        dens_x = create_r(setfl.n_r*setfl.d_r,setfl.n_r)
        dens_y = setfl.density_function(s)

        ax[1].plot(dens_x,dens_y,label=s)
    ax[1].legend()

    # plot the embedding energy function
    for s in setfl.elements:
        embed_x = create_r(setfl.n_rho*setfl.d_rho,setfl.n_rho)
        embed_y = setfl.embedding_function(s)
        ax[2].plot(embed_x,embed_y,label=s)
    ax[2].legend()

    fig.tight_layout()

    plt.show()
if __name__ == "__main__":

    #setfl_filename = os.path.join(
    #        'test_EamSetflFile',
    #        'Ni1_Mendelev_2010.eam.fs')

    setfl_filename = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')

    setfl_reader = SeatonSetflReader(path=setfl_filename)
    setfl_reader.read()
    print("elements: ", setfl_reader.elements)
    print("element_pairs: ", setfl_reader.element_pairs)
    print("n_rho: ", setfl_reader.n_rho)
    print("d_rho: ", setfl_reader.d_rho)
    print("n_r: ", setfl_reader.n_r)
    print("d_r: ", setfl_reader.d_r)
    print("cutoff: ", setfl_reader.cutoff)
    assert len(setfl_reader.embedding_function("Ni")) == setfl_reader.n_r
    assert len(setfl_reader.density_function("Ni")) == setfl_reader.n_rho
    assert len(setfl_reader.pair_function("NiNi")) == setfl_reader.n_r

    plot_eam_potential_from_setflfile(setfl_filename)
