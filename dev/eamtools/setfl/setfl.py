#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = "Eugene Ragasa"
__copyright__ = "Copyright 2019, Eugene J. Ragasa"
__licence__ = "MIT"
__contributors__ = ["Eugene J Ragasa, R. Seaton Ullberg, Simon Phillpot"]
__maintainer__ = "Eugene J. Ragasa"

def get_setfl_pair_order(symbols):
    pairs = []
    for i,s_i in enumerate(symbols):
        for j,s_j in enumerate(symbols):
            if i >= j:
                pairs.append([s_i,s_j])
    return pairs
