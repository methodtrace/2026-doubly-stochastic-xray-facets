"""Families of facet rays and fan rays: facet fibers, the quadratic family, refinement descent, two-point block lifts,
padded cubical corners and the binary interval-cancellation certificate. Exact (integers and Fractions); graphs and
ranks from exact_graphs.py (forced-edge matching, rational elimination).

Potentials u on {0..k-1}; v(t) = min_{i+j=t} u(i) + u(j); E(u) the tight graph; A = E*(u) the allowed graph;
Z(G), Y(G), bip(G) as in facet_criterion.py; spr_m u(i) = u(ceil(i/m)) on {0..m(k-1)}.

  cubical_fiber(m)         chi_m(i) = chi(floor(i/m)) on k = 19m points (chi the k = 19 potential of ex:noTS):
                           v = 1 exactly on {6m-1, 11m-1, 14m-1, 28m-1}; 19 m^2 tight edges are not allowed; dim Z = 3,
                           dim Z(A) = 3 + m, bip(A) = m, face dimension 2k - 4 (a facet without total support); the
                           bipartite components are the pairs I_a = 11m + a, J_a = 17m - 1 - a (a < m), J_a the only
                           tight neighbour of I_a; the fiber {u + sum_a lambda_a (e_{I_a} - e_{J_a}) : lambda in
                           [-1, 0]^m} consists of potentials defining the same inequality; its 2^m vertices are rays
                           (dim Z = 3) and its centre has tight graph A (total support).
  quadratic(a, b)          G_ab with increments (b^a, 0^(b-1), (b-1)a + 1, 1^b), n = a + 2b (k = n + 1): the stated v
                           and tight graph, a tight permutation fixing e = a + b - 1 (universal vertex), total support,
                           connected; dim Z = 3 iff a <= b + 1, else 4; height (b - 1)a + 1; P_L, P_R; core.
  descent(u, m)            sparse refinement: refined tight edges project to base tight edges and (ma, mb) lifts
                           (a, b); total support and connectivity descend, with a block-collapse doubly stochastic
                           certificate D = C/m + (1 - 1/m) e_0 e_0^T supported on E(u) for every base edge; under the
                           flat-level common-difference condition, ray descends.
  two_point_lifts(w)       the extensions u = (-p - q, -q, w) of a normalized facet ray w keeping its tight pairs:
                           p + q <= A, q <= B; rays exactly on two of the lines p = q, p = s0, p + q = A, q = B.
  padded_corners(m, a, b)  pad_L^a pad_R^b of the normalized corners of the cubical fiber: height-two rays.
  binary_cancellation(d)   the interval-cancellation closure on a binary increment word; two final classes certify a
                           ray (sound; complete for k <= 12 by exhaustion).
"""
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd

from . import exact_graphs as GR

CHI_ZEROS_COMPLEMENT = {2, 4, 5, 7, 9, 10, 11, 12, 13}


# ---------------- facet fibers: cubical family ----------------
def chi(m=1):
    return [int(i // m in CHI_ZEROS_COMPLEMENT) for i in range(19 * m)]


def cubical_fiber(m):
    """Checks all statements listed in the module docstring for chi_m; returns a summary dict (raises
    AssertionError on failure)."""
    u = chi(m); k = len(u)
    v, tight = GR.tight_graph(u)
    assert v == [int(t in {6 * m - 1, 11 * m - 1, 14 * m - 1, 28 * m - 1}) for t in range(2 * k - 1)]
    A = GR.allowed_edges(k, tight)
    assert A is not None and len(tight - A) == 19 * m * m
    d = GR.dimensions(k, A); dt = GR.dimensions(k, tight, face=False)
    assert dt["zdim"] == 3 and d["zdim"] == 3 + m and d["ydim"] == m and d["face_dim"] == 2 * k - 4
    comps = GR.components(k, A)
    pairs = [(11 * m + a, 17 * m - 1 - a) for a in range(m)]
    assert sorted(c["vertices"] for c in comps if c["bipartite"]) == sorted(sorted(p) for p in pairs)
    assert {i + j for i, j in A} == set(range(2 * k - 1))
    for i, j in pairs:
        assert sorted(b for a, b in tight if a == i) == [j] and v[2 * i] == v[2 * j] == 0
    rays = 0
    for lam in product([-1, 0], repeat=m):
        w = u[:]
        for (i, j), x in zip(pairs, lam):
            w[i] += x; w[j] -= x
        assert sum(w) == sum(u)
        vw, Ew = GR.tight_graph(w)
        assert vw == v and A <= Ew and GR.zdim(k, Ew) == 3
        rays += 1
    c = [Fraction(x) for x in u]
    for i, j in pairs:
        c[i] -= Fraction(1, 2); c[j] += Fraction(1, 2)
    vc, Ec = GR.tight_graph(c)
    assert vc == v and Ec == A
    return {"m": m, "k": k, "lost_edges": len(tight - A), "zdim": 3, "zdim_allowed": d["zdim"], "bip": d["ydim"],
            "face_dim": d["face_dim"], "fiber_vertices_rays": rays}


def zero_set_certificate():
    """The zero set S = {0,1,3,6,8,14,15,17,18} of chi has S + S = [0, 36] minus {5, 10, 13, 27}."""
    S = set(range(19)) - CHI_ZEROS_COMPLEMENT
    return {a + b for a in S for b in S} == set(range(37)) - {5, 10, 13, 27}


# ---------------- the quadratic family ----------------
def quadratic_potential(a, b):
    d = [b] * a + [0] * (b - 1) + [(b - 1) * a + 1] + [1] * b
    u = GR.from_increments(d)
    return d, u, [x - i for i, x in enumerate(u)]


def quadratic_claimed_v(a, b):
    p, e, n = b - 1, a + b - 1, a + 2 * b
    L, C, T, S = p * (a - 1), p * (2 * a - 1), a + e, e + 1 + n
    out = []
    for t in range(2 * n + 1):
        if t <= a:
            x = p * t
        elif t <= e:
            x = a * b - t
        elif t <= T:
            x = L + p * (t - e)
        elif t <= 2 * e:
            x = 2 * a * b - t
        elif t < S:
            x = L + C
        else:
            x = 2 * C
        out.append(x)
    return out


def quadratic_claimed_edges(a, b):
    e, n = a + b - 1, a + 2 * b
    A, B, Q = range(a + 1), range(a, e + 1), range(e + 1, n + 1)
    T, S = a + e, e + 1 + n
    E = {(i, j) for i in A for j in A if i + j <= a}
    E |= {(0, j) for j in B} | {(i, e) for i in A}
    E |= {(i, j) for i in B for j in B if i + j >= T}
    E |= {(e, j) for j in Q} | {(i, j) for i in Q for j in Q if i + j >= S}
    E.add((a - 1, n))
    if a <= b + 1:
        E.add((0, T))
    return E | {(j, i) for i, j in E}


def quadratic_permutation(a, b):
    e, n = a + b - 1, a + 2 * b
    T, S = a + e, e + 1 + n
    return [a - i if i <= a else T - i if i < e else e if i == e else S - i for i in range(n + 1)]


def quadratic(a, b, forced_edges=False):
    """Checks the theorem for G_ab; returns (dim Z, height, k). forced_edges: also every tight edge allowed by
    forced matching, and the face dimension 2k - 1 - dim Z by rank."""
    d, u, w = quadratic_potential(a, b)
    k, e, n = len(w), a + b - 1, a + 2 * b
    v, E = GR.tight_graph(w)
    assert v == quadratic_claimed_v(a, b) and E == quadratic_claimed_edges(a, b)
    s = quadratic_permutation(a, b)
    assert sorted(s) == list(range(k)) and s[e] == e and all((i, s[i]) in E for i in range(k))
    assert all((e, j) in E for j in range(k)) and GR.bipartite_connected(k, E)
    zd = GR.zdim(k, E)
    assert zd == (3 if a <= b + 1 else 4)
    H = (b - 1) * a + 1
    assert min(d) == 0 and reduce(gcd, d) == 1 and max(d) == H
    assert d[a] == 0 and (0, a) in E and d[e] == H and (e + 1, n) in E and 0 < d[0] and d[-1] < H
    if forced_edges:
        assert GR.allowed_edges(k, E) == E
        assert GR.dimensions(k, E)["face_dim"] == 2 * k - 1 - zd
    return zd, H, k


# ---------------- sparse refinement ----------------
def spr(u, m):
    return [u[(i + m - 1) // m] for i in range(m * (len(u) - 1) + 1)]


def total_support(k, E):
    A = GR.allowed_edges(k, E)
    return A is not None and A == E


def lifting_condition(u, v, E):
    """The flat-level common-difference condition, as a rank condition on Z(u)."""
    k = len(u)
    canon = []
    for s in range(2 * k - 1):
        pr = [(i, s - i) for i in range(k) if i <= s - i < k and (i, s - i) in E]
        row = [0] * k; i, j = pr[0]; row[i] += 1; row[j] += 1
        canon.append(row)
    rows = GR.z_rows(k, E)
    flat = [s for s in range(2, 2 * k - 1) if v[s] == v[s - 1] and any(i > 0 and j > 0 and i + j == s for i, j in E)]
    diffs = [[x - y for x, y in zip(canon[s], canon[s - 1])] for s in flat]
    eqs = [[x - y for x, y in zip(dd, diffs[0])] for dd in diffs[1:]]
    return GR.rank(rows + eqs, k) == GR.rank(rows, k)


def collapse_certificates(u, m, E, Ep):
    """For every base tight edge (a, b): a refined tight permutation through (ma, mb) and the matrix
    m D = C + (m - 1) e_0 e_0^T (C the block counts), checked doubly stochastic (row and column sums m), supported
    on E(u) and positive at (a, b). Returns the number of certificates."""
    k = len(u); kp = m * (k - 1) + 1; n = 0
    for a, b in sorted(E):
        pi = GR.matching(kp, Ep, (m * a, m * b))
        assert pi is not None and pi[m * a] == m * b
        N = [[0] * k for _ in range(k)]
        for i, j in enumerate(pi):
            N[0 if i == 0 else (i - 1) // m + 1][0 if j == 0 else (j - 1) // m + 1] += 1
        N[0][0] += m - 1
        assert all(sum(r) == m for r in N) and all(sum(r[j] for r in N) == m for j in range(k))
        assert N[a][b] > 0 and all(not N[i][j] or (i, j) in E for i in range(k) for j in range(k))
        n += 1
    return n


def descent(u, m):
    """Returns (base total support, refined total support, base connected, refined connected, base dim Z,
    refined dim Z, lifting condition, number of collapse certificates); asserts the descent statements."""
    k = len(u); up = spr(u, m); kp = len(up)
    v, E = GR.tight_graph(u); vp, Ep = GR.tight_graph(up)
    assert all(((i + m - 1) // m, (j + m - 1) // m) in E for i, j in Ep)
    assert all((m * a, m * b) in Ep for a, b in E)
    bts, rts = total_support(k, E), total_support(kp, Ep)
    bc, rc = GR.bipartite_connected(k, E), GR.bipartite_connected(kp, Ep)
    assert (not rc or bc) and (not rts or bts)
    ncert = collapse_certificates(u, m, E, Ep) if rts else 0
    bz, rz = GR.zdim(k, E), GR.zdim(kp, Ep)
    cond = lifting_condition(u, v, E)
    assert not (cond and rz == 3) or bz == 3
    return bts, rts, bc, rc, bz, rz, cond, ncert


# ---------------- ray-only padding ----------------
def padding_witnesses(u, E):
    d = [y - x for x, y in zip(u, u[1:])]; h = max(d)
    pl = [q for q in range(len(d)) if d[q] == 0 and (0, q) in E]
    pr = [q for q in range(len(d)) if d[q] == h and (q + 1, len(u) - 1) in E]
    return pl, pr


def padding_ray_only(word):
    """For a normalized ray with increments word: pad_L is a ray iff P_L, pad_R is a ray iff P_R. Returns
    (is ray, P_L, P_R) and asserts the equivalences."""
    u = GR.from_increments(word); k = len(u)
    _, E = GR.tight_graph(u)
    if GR.zdim(k, E) != 3:
        return False, None, None
    pl, pr = padding_witnesses(u, E)
    ul, ur = [u[0]] + u, u + [u[-1] + max(word)]
    assert (GR.zdim(k + 1, GR.tight_graph(ul)[1]) == 3) == bool(pl)
    assert (GR.zdim(k + 1, GR.tight_graph(ur)[1]) == 3) == bool(pr)
    return True, bool(pl), bool(pr)


# ---------------- two-point block lifts ----------------
def pattern(u):
    d = [Fraction(b) - Fraction(a) for a, b in zip(u, u[1:])]
    lo = min(d); d = [x - lo for x in d]
    den = 1
    for x in d:
        den = den * x.denominator // gcd(den, x.denominator)
    d = [int(x * den) for x in d]; g = reduce(gcd, d)
    return tuple(x // g for x in d)


def lift_candidates(w):
    """(A, B, s0, points): the at most six intersection points of two of the lines with p + q <= A, q <= B."""
    n = len(w); v = GR.tight_graph(w)[0]; s0 = w[1]
    A = min(w[t + 2] - v[t] for t in range(n - 2)); B = min(w[t + 1] - v[t] for t in range(n - 1))
    pts = {(s0, s0), (Fraction(A, 2), Fraction(A, 2)), (B, B), (s0, A - s0), (s0, B), (A - B, B)}
    return A, B, s0, sorted((Fraction(p), Fraction(q)) for p, q in pts if p + q <= A and q <= B)


def ray_ts_connected(u):
    k = len(u); _, E = GR.tight_graph(u)
    return GR.zdim(k, E) == 3, total_support(k, E), GR.bipartite_connected(k, E), E


def two_point_lifts(w, grid=False):
    """Children patterns of the normalized facet ray w (as a list of values with w[0] = 0); asserts that each child
    is a ray with total support and connected tight graph keeping the shifted tight pairs of w. grid=True also tests
    every point of a half-integer grid: (L) iff the shifted pairs are kept, and ray iff on at least two lines.
    Returns (children, grid points tested)."""
    A, B, s0, pts = lift_candidates(w)
    _, Ew = GR.tight_graph(w)
    out = []
    for p, q in pts:
        u = [-p - q, -q] + list(w)
        ray, ts, conn, E = ray_ts_connected(u)
        assert ray and ts and conn and all((i + 2, j + 2) in E for i, j in Ew)
        out.append(pattern(u))
    tested = 0
    if grid:
        for p2 in range(-4, 2 * int(max(A, s0)) + 6):
            for q2 in range(-4, 2 * int(max(B, s0)) + 6):
                p, q = Fraction(p2, 2), Fraction(q2, 2)
                u = [-p - q, -q] + list(w)
                _, E = GR.tight_graph(u)
                kept = all((i + 2, j + 2) in E for i, j in Ew)
                assert kept == (p + q <= A and q <= B)
                if kept:
                    on = (p == q) + (p == s0) + (p + q == A) + (q == B)
                    assert (GR.zdim(len(u), E) == 3) == (on >= 2)
                tested += 1
    return out, tested


def mirror(p):
    H = max(p)
    return tuple(H - x for x in reversed(p))


# ---------------- padded cubical corners (height-two rays) ----------------
def cubical_corners(m):
    u = chi(m); out = []
    pairs = [(11 * m + a, 17 * m - 1 - a) for a in range(m)]
    for lam in product([0, -1], repeat=m):
        f = u[:]
        for (i, j), x in zip(pairs, lam):
            f[i] += x; f[j] -= x
        out.append(pattern([f[i] + i for i in range(len(f))]))
    return out


def padded_corner_words(k):
    """All pad_L^a pad_R^b of normalized corners with 19m + a + b = k, m >= 1."""
    words = []
    for m in range(1, k // 19 + 1):
        for d in cubical_corners(m):
            for a in range(k - 19 * m + 1):
                b = k - 19 * m - a
                words.append((m, a, b, (0,) * a + tuple(d) + (max(d),) * b))
    return words


# ---------------- binary interval cancellation ----------------
def binary_cancellation(d):
    """(dim Z, closure gives two classes): for increments d in {0, 1}; the closure merges increment positions using
    the balance identities of consecutive tight pairs on a level (each a signed sum of increments with zero sum)
    whenever exactly two classes survive with nonzero coefficients."""
    n = len(d); k = n + 1; u = GR.from_increments(d)
    rows = []
    for t in range(2 * k - 1):
        prs = [(p, t - p) for p in range(max(0, t - k + 1), t // 2 + 1)]
        mn = min(u[p] + u[q] for p, q in prs)
        tight = [(p, q) for p, q in prs if u[p] + u[q] == mn]
        for (a, _), (b, _) in zip(tight, tight[1:]):
            row = [0] * n
            for i in range(a, b):
                row[i] += 1
            for i in range(t - b, t - a):
                row[i] -= 1
            assert sum(row) == 0 and sum(x * y for x, y in zip(row, d)) == 0
            rows.append(row)
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    changed = True
    while changed:
        changed = False
        for row in rows:
            co = {}
            for i, x in enumerate(row):
                if x:
                    c = find(i); co[c] = co.get(c, 0) + x
            co = {c: x for c, x in co.items() if x}
            if len(co) == 2:
                a, b = co
                assert co[a] == -co[b] and d[a] == d[b]
                par[find(a)] = find(b); changed = True
    classes = len({find(i) for i in range(n)})
    zd = 1 + n - GR.rank(rows, n)
    return zd, len(set(d)) > 1 and classes == 2
