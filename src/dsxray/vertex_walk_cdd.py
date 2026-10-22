"""Independent recomputation of the fan rays: walk on the vertices of Q = sum_t Q_t with cddlib normal cones.

Requires the optional dependency pycddlib (import name `cdd`, version 3 API). Shares no code with fan_walk.py: no
difference matroid, no flat enumeration, no ridge crossing.

A vertex of Q is a choice of one pair (i, t-i), i <= t-i, on every anti-diagonal 1 <= t <= 2k-3 that is the unique
minimizer of w(i) + w(j) for some w. Its normal cone, in coordinates w(0) = w(k-1) = 0, is {w : w(a)+w(b) <= w(c)+w(d)}
for the chosen (a, b) and every other pair (c, d) of the same anti-diagonal. cddlib (floating point) gives the extreme
rays; each is rounded to a primitive integer vector and verified exactly: it satisfies every inequality of its cone,
and its tight inequalities have rank k-3. The walls (facets of the cone) are the rows whose tight rays have rank k-3
(floating-point rank of a small integer matrix);
crossing the wall of anti-diagonal t replaces the chosen pair of t. Every facet of Q contains a vertex and the vertex
graph is connected, so the union of the rays of all vertex cones is the set of all rays. A ray or wall missed by the
floating-point computation would go unnoticed; every reported ray is correct.
"""
import math
import random
import time
from fractions import Fraction
from functools import reduce
from math import gcd

import numpy as np


def _nint(x):
    """Nearest integer (halves upward)."""
    return math.floor(x + 0.5)


def _cdd():
    import cdd
    return cdd


def pairs(k, t):
    return [(i, t - i) for i in range(k) if i <= t - i < k]


def ineqs(k, choice):
    rows = []; ids = []
    for t in range(1, 2 * k - 2):
        a, b = choice[t]
        for (c, d) in pairs(k, t):
            if (c, d) == (a, b):
                continue
            v = [0] * k; v[c] += 1; v[d] += 1; v[a] -= 1; v[b] -= 1
            rows.append(v[1:k - 1]); ids.append((t, (c, d)))
    return rows, ids


def iprim(v):
    g = reduce(gcd, [abs(x) for x in v]) or 1
    return tuple(x // g for x in v)


def rank_exact(rows):
    M = [[Fraction(x) for x in r] for r in rows]; rk = 0; n = len(M[0]) if M else 0
    for c in range(n):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def cone_data(k, choice, verified):
    """(walls, rays) of the normal cone of a vertex; every new ray is verified exactly and added to `verified`."""
    cdd = _cdd()
    rows, ids = ineqs(k, choice)
    P = cdd.polyhedron_from_matrix(cdd.matrix_from_array([[0] + r for r in rows], rep_type=cdd.RepType.INEQUALITY))
    gen = cdd.copy_generators(P).array
    rays = []
    for g in gen:
        if abs(g[0]) > 1e-9:
            continue                                  # the apex
        vec = list(g[1:]); m = max(abs(x) for x in vec); vec = [x / m for x in vec]
        for mult in range(1, 100000):
            iv = [_nint(x * mult) for x in vec]
            if all(abs(x * mult - y) < 1e-7 for x, y in zip(vec, iv)):
                break
        else:
            raise RuntimeError("could not rationalize a ray")
        rays.append(iprim(iv))
    R = np.array(rays, dtype=np.int64); A = np.array(rows, dtype=np.int64)
    V = A @ R.T                                       # exact integer incidence
    if V.min() < 0:
        raise RuntimeError("ray violates its cone")
    walls = []
    for i in range(len(rows)):
        T = R[V[i] == 0]
        if len(T) >= k - 3 and np.linalg.matrix_rank(T.astype(float)) == k - 3:
            walls.append(ids[i])
    for j, r in enumerate(rays):
        if r in verified:
            continue
        if rank_exact([rows[i] for i in range(len(rows)) if V[i, j] == 0]) != k - 3:
            raise RuntimeError("reported ray is not extreme")
        verified.add(r)
    return walls, rays


def vertex_of(k, w):
    ch = {}
    for t in range(1, 2 * k - 2):
        vals = sorted((w[i] + w[j], (i, j)) for (i, j) in pairs(k, t))
        if len(vals) > 1 and vals[0][0] == vals[1][0]:
            return None
        ch[t] = vals[0][1]
    return ch


_K = None
_VERIFIED = set()


def _init(k):
    global _K, _VERIFIED
    _K = k; _VERIFIED = set()


def _expand(keys):
    k = _K; T = range(1, 2 * k - 2); out = []
    for b in keys:
        ch = {t: (i, t - i) for t, i in zip(T, b)}
        walls, rays = cone_data(k, ch, _VERIFIED)
        nbs = []
        for (t, cd) in walls:
            nb = dict(ch); nb[t] = cd; nbs.append(bytes(nb[s][0] for s in T))
        out.append(([iprim((0,) + r + (0,)) for r in rays], nbs))
    return out


def walk(k, workers=1, seed=1, log=None, chunk=256):
    """Level-synchronous walk of the vertex graph of Q. Returns (number of vertices, set of rays)."""
    rng = random.Random(seed)
    while True:
        w = [0] + [rng.randint(-10**9, 10**9) for _ in range(k - 2)] + [0]
        ch = vertex_of(k, w)
        if ch:
            break
    b0 = bytes(ch[t][0] for t in range(1, 2 * k - 2))
    seen = {b0}; frontier = [b0]; rays = set(); nv = 0; level = 0; t0 = time.time()
    pool = None
    if workers > 1:
        from multiprocessing import Pool
        pool = Pool(workers, initializer=_init, initargs=(k,))
    else:
        _init(k)
    try:
        while frontier:
            chunks = [frontier[i:i + chunk] for i in range(0, len(frontier), chunk)]
            results = pool.imap_unordered(_expand, chunks) if pool else map(_expand, chunks)
            nxt = []
            for res in results:
                for rs, nbs in res:
                    nv += 1; rays.update(rs)
                    for kk in nbs:
                        if kk not in seen:
                            seen.add(kk); nxt.append(kk)
            frontier = nxt; level += 1
            if log:
                print(f"  k={k} level {level} vertices {nv} next {len(frontier)} rays {len(rays)} "
                      f"{time.time() - t0:.0f}s", file=log, flush=True)
    finally:
        if pool:
            pool.close(); pool.join()
    return nv, rays
