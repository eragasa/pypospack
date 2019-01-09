# -*- coding: utf-8 -*-
"""This module provides pareto calculation functions

Eugene J. Ragasa, University of Florida, developed the original version
of this code.  Dmitriy Morozov, Lawrence Berkeley Labs, provided speed ups for the
Pareto versions of the code in Dec 2016.


"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"
import os
#import pareto

from random     import randint, seed
import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('%s function took %0.3f ms' % (f.__name__, (time2-time1)*1000.0))
        return ret
    return wrap

def dominates(p1, p2):
    for x,y in zip(p1,p2):
        if y < x:
            return False
    return True

def pareto_bruteforce(pts, indices = None):
    if indices is None:
        indices = list(range(len(pts)))
    result = []
    for i in indices:
        for j in indices:
            if i == j: continue
            if dominates(pts[j], pts[i]):
                break
        else:
            result.append(i)
    return result

def pareto_merge(lo, hi, i, dim):
    if len(lo) == 0 or len(hi) == 0:
        return lo + hi

    survivors = set()
    for j in range(dim):
        if i == j: continue
        m = min(p[j] for p in lo)
        survivors.update(k for k in range(len(hi)) if hi[k][j] < m)

    return lo + [hi[k] for k in survivors]

def pareto(pts,chunk_sz=500,is_debug=False):
    idx = list(range(len(pts)))
    sz = chunk_sz

    old_len_idx = len(idx)
    while True:
        if len(idx) < sz:
            break
        if is_debug: print('total_idx {} > chunk_sz {}'.format(len(idx),sz))
        chunks=zip(*[iter(idx)]*sz)
        idx= []
        for i,c in enumerate(chunks):
            #print("\tworking on chunk {}".format(i))
            idx += pareto_bruteforce(pts,c)

        new_len_idx = len(idx)
        if (old_len_idx == new_len_idx):
            break
        if (old_len_idx - new_len_idx) < sz:
            sz = 2*sz
        old_len_idx = new_len_idx
    if is_debug: print("running final pareto calculation on {} points.".format(len(idx)))
    return pareto_bruteforce(pts,idx)

def pareto_lawrenceberkeley(pts, indices = None, i = 0):
    if indices is None:
        indices = list(range(len(pts)))
    l = len(indices)
    if l <= 1:
        return indices

    if l < 1000:
        return pareto_bruteforce(pts, indices)

    indices.sort(key = lambda x: pts[x][i])     # lazy: should use partition instead

    dim = len(pts[0])
    optimalLo = pareto(pts, indices[:l//2], (i + 1) % dim)
    optimalHi = pareto(pts, indices[l//2:], (i + 1) % dim)

    return pareto_bruteforce(pts, optimalLo + optimalHi)     # lazy: FIXME
    #return pareto_merge(optimalLo, optimalHi, i, dim)

def read_data(fn):
    with open(fn) as f:
        lines = f.readlines()
        names = ['sim_id'] + [x.split() for x in lines[0].split('|')]
        names[1] = names[1][1:]     # skip sim_id
        values = []
        for line in lines[1:]:
            line = line.split('|')
            line = [x.split() for x in line]
            values.append([
                int(line[0][0]),
                [float(x) for x in line[0][1:]],
                [float(x) for x in line[1]]])

        return names, values
