"""Independent completeness check for the fan rays, exact version (no floating point, no dependencies).

Walks the vertex graph of Q = sum_t Q_t. A vertex is a choice of one pair (i, t-i), i <= t-i, on every anti-diagonal
1 <= t <= 2k-3 that is the unique minimizer of w(i)+w(j) for some w. Its normal cone, in coordinates w(0)=w(k-1)=0
(which removes the lineality space spanned by 1 and iota), is {w in Z^(k-2) : a . w >= 0 for every row a}, one row per
non-chosen pair (c,d) on each anti-diagonal: a = e_c + e_d - e_a - e_b restricted to coordinates 1..k-2. The extreme rays
of every vertex cone are computed by the double description method in exact integer arithmetic; the walls of the cone
are the rows whose tight rays have rank k-3; crossing the wall of anti-diagonal t to the pair (c,d) gives the neighboring
vertex. Every facet of Q contains a vertex and the vertex graph of a polytope is connected (Balinski), so the union of the
rays of all vertex cones is the set of all rays of the fan. This module shares no code with fan_walk.py.
"""
import random
import time
from math import gcd
from functools import reduce


def pairs(k, t):
    return [(i, t - i) for i in range(k) if i <= t - i < k]


def iprim(v):
    g = reduce(gcd, (abs(x) for x in v), 0) or 1
    return tuple(x // g for x in v)


def rank_int(rows):
    """Exact rank of an integer matrix (fraction-free elimination on Python ints)."""
    M = [list(r) for r in rows]
    if not M: return 0
    n = len(M[0]); rk = 0
    for c in range(n):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; piv = M[rk]
        for i in range(rk + 1, len(M)):
            if M[i][c]:
                f, g = M[i][c], piv[c]
                row = [g * x - f * y for x, y in zip(M[i], piv)]
                h = reduce(gcd, (abs(x) for x in row), 0) or 1
                M[i] = [x // h for x in row]
        rk += 1
    return rk


def cone_rows(k, choice):
    rows, ids = [], []
    for t in range(1, 2 * k - 2):
        a, b = choice[t]
        for (c, d) in pairs(k, t):
            if (c, d) == (a, b): continue
            v = [0] * k; v[c] += 1; v[d] += 1; v[a] -= 1; v[b] -= 1
            rows.append(tuple(v[1:k - 1])); ids.append((t, (c, d)))
    return rows, ids


def dot(a, b): return sum(x * y for x, y in zip(a, b))


def sparse(a): return tuple((j, x) for j, x in enumerate(a) if x)


def sdot(sa, r): return sum(x * r[j] for j, x in sa)


def solve_basis(B):
    """Integer vectors r_j with B r_j = c_j e_j, c_j > 0 (columns of a positive multiple of B^{-1}); B square nonsingular.
    Integer Gauss-Jordan on [B | I]: rows are only combined and divided by row gcds, so every row stays a nonzero
    multiple of a row of the rational reduced form."""
    n = len(B)
    M = [list(B[i]) + [int(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]; pv = M[c][c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                row = [pv * x - f * y for x, y in zip(M[i], M[c])]
                g = reduce(gcd, (abs(x) for x in row), 0) or 1
                M[i] = [x // g for x in row]
    # now M = [D | X] with D diagonal, and B^{-1} row i = X_i / D_ii
    D = [M[i][i] for i in range(n)]
    L = reduce(lambda a, b: a * b // gcd(a, b), (abs(d) for d in D), 1)
    out = []
    for j in range(n):
        col = [M[i][n + j] * (L // D[i]) for i in range(n)]      # L * (B^{-1})_{ij}, exact since D[i] | L (up to sign)
        out.append(iprim(col))
    return out


def independent_rows(rows, d):
    """Indices of the first d linearly independent rows (greedy), by incremental integer elimination."""
    basis = []; piv = []; idx = []
    for i, r in enumerate(rows):
        v = list(r)
        for (c, b) in zip(piv, basis):
            if v[c]:
                f, g = v[c], b[c]; v = [g * x - f * y for x, y in zip(v, b)]
                h = reduce(gcd, (abs(x) for x in v), 0) or 1; v = [x // h for x in v]
        c = next((j for j, x in enumerate(v) if x), None)
        if c is None: continue
        basis.append(v); piv.append(c); idx.append(i)
        if len(idx) == d: break
    return idx


def extreme_rays(rows, d):
    """Double description for the pointed full-dimensional cone {w in R^d : row . w >= 0}. Returns primitive integer rays."""
    # initial basis: greedy independent rows
    bidx = independent_rows(rows, d)
    if len(bidx) < d: raise RuntimeError("cone is not pointed")
    R = solve_basis([rows[i] for i in bidx])                         # B r_j = c e_j with c > 0: row j positive on r_j, others zero
    rays = [r for r in R]
    done = list(bidx)
    sp = [sparse(a) for a in rows]
    zero = [frozenset(i for i in bidx if sdot(sp[i], r) == 0) for r in rays]
    for i, a in enumerate(rows):
        if i in bidx: continue
        sa = sp[i]
        vals = [sdot(sa, r) for r in rays]
        P = [j for j, v in enumerate(vals) if v > 0]; N = [j for j, v in enumerate(vals) if v < 0]; Z = [j for j, v in enumerate(vals) if v == 0]
        new, newz = [], []
        if N:
            for p in P:
                for q in N:
                    common = zero[p] & zero[q]
                    if len(common) < d - 2: continue
                    # combinatorial adjacency: no third ray is tight on all of `common`
                    if any(common <= zero[s] for s in range(len(rays)) if s != p and s != q): continue
                    r = iprim([vals[p] * y - vals[q] * x for x, y in zip(rays[p], rays[q])])
                    new.append(r); newz.append(common | {i})
        keep = P + Z
        rays = [rays[j] for j in keep] + new
        zero = [zero[j] | ({i} if vals[j] == 0 else set()) for j in keep] + newz
        zero = [frozenset(z) for z in zero]
        done.append(i)
    return rays


def cone_data(k, choice):
    rows, ids = cone_rows(k, choice)
    d = k - 2
    rays = extreme_rays(rows, d)
    walls = []
    for a, idt in zip(rows, ids):
        sa = sparse(a)
        T = [r for r in rays if sdot(sa, r) == 0]
        if len(T) >= d - 1 and rank_int(T) == d - 1: walls.append(idt)
    return walls, rays


def vertex_of(k, w):
    ch = {}
    for t in range(1, 2 * k - 2):
        vals = sorted((w[i] + w[j], (i, j)) for (i, j) in pairs(k, t))
        if len(vals) > 1 and vals[0][0] == vals[1][0]: return None
        ch[t] = vals[0][1]
    return ch


_K = None
def _init(k):
    global _K; _K = k


def _expand(keys):
    k = _K; T = range(1, 2 * k - 2); out = []
    for b in keys:
        ch = {t: (i, t - i) for t, i in zip(T, b)}
        walls, rays = cone_data(k, ch)
        nbs = []
        for (t, cd) in walls:
            nb = dict(ch); nb[t] = cd; nbs.append(bytes(nb[s][0] for s in T))
        out.append(([(0,) + r + (0,) for r in rays], nbs))
    return out


def walk(k, workers=1, seed=1, log=None, chunk=64):
    """Returns (number of vertices of Q, set of rays normalized by w(0) = w(k-1) = 0)."""
    rng = random.Random(seed)
    while True:
        w = [0] + [rng.randint(-10**9, 10**9) for _ in range(k - 2)] + [0]
        ch = vertex_of(k, w)
        if ch: break
    T = range(1, 2 * k - 2)
    b0 = bytes(ch[t][0] for t in T)
    seen = {b0}; frontier = [b0]; rays = set(); nv = 0; t0 = time.time(); level = 0
    pool = None
    if workers > 1:
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(k,))
    else:
        _init(k)
    try:
        while frontier:
            chunks = [frontier[i:i + chunk] for i in range(0, len(frontier), chunk)]
            results = pool.imap_unordered(_expand, chunks) if pool else map(_expand, chunks)
            nxt = []
            for res in results:
                for rs, nbs in res:
                    nv += 1
                    for r in rs: rays.add(iprim(r))
                    for kk in nbs:
                        if kk not in seen: seen.add(kk); nxt.append(kk)
            level += 1
            if log:
                print(f"  k={k} level {level} vertices {nv} next {len(nxt)} rays {len(rays)} {time.time()-t0:.0f}s", file=log, flush=True)
            frontier = nxt
    finally:
        if pool: pool.close(); pool.join()
    return nv, rays
