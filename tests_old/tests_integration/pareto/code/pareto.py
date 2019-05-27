from __future__ import print_function
from builtins   import zip, range

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

def pareto(pts, indices = None, i = 0):
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


if __name__ == '__main__':
    seed(10)

    pts = []
    for i in range(5000):
        pts.append(tuple(randint(0,100) for _ in range(8)))

    res = timing(pareto_bruteforce)(pts)
    print(len(res))

    res2 = timing(pareto)(pts)
    print(len(res2))
