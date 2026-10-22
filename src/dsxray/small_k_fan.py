"""Independent enumerations for small k: line enumeration of the fan and exact facet lists from the profiles.

Line enumeration. The fan is cut out by the comparison hyperplanes u(i) + u(j) = u(i') + u(j'), i + j = i' + j'
(coordinates u(1..k-2), u(0) = u(k-1) = 0). Every ray lies on a line that is the intersection of k-3 independent
comparison hyperplanes; such a line is the kernel of the (k-3) x (k-2) matrix of their normals, computed as the
vector of signed maximal minors (integer exterior product). A direction d (either sign) is a ray iff the minimizing
pairs of every anti-diagonal of d span a space of rank k-3 (exact rank). This adapts the edge-arrangement method of
Gritzmann and Sturmfels and shares no code with fan_walk.py.

Profiles. hat-P(k) is the convex hull of the profiles x(sigma)_t = #{i : i + sigma(i) = t} of all permutations. For
k <= 8 every candidate inequality (from a ray or a coordinate) is checked against all profiles in integer arithmetic:
valid on all, and tight on an affinely independent set of 2k - 3 profiles (exact rank). This gives exact facet lists
whose completeness rests on the fact that every facet other than a coordinate one comes from a ray.
"""
import itertools
from fractions import Fraction as F
from functools import reduce
from math import gcd

import numpy as np


def primitive(v):
    vals = [F(x) for x in v]
    den = reduce(lambda a, b: a * b // gcd(a, b), (x.denominator for x in vals), 1)
    out = [int(x * den) for x in vals]; g = reduce(gcd, (abs(x) for x in out), 0)
    return tuple(x // g for x in out) if g else tuple(out)


def int_rank(rows, stop=None):
    """Exact rank of an integer matrix by fraction-free elimination; stops early once `stop` is reached."""
    basis = []
    for r in rows:
        v = [int(x) for x in r]
        for (c, b) in basis:
            if v[c]:
                f, gg = v[c], b[c]
                v = [gg * x - f * y for x, y in zip(v, b)]
                h = reduce(gcd, (abs(x) for x in v), 0) or 1
                v = [x // h for x in v]
        c = next((j for j, x in enumerate(v) if x), None)
        if c is None:
            continue
        basis.append((c, v))
        if stop is not None and len(basis) >= stop:
            break
    return len(basis)


# ---------------- line enumeration ----------------
def fan_setup(k):
    d = k - 2; groups = []; comparisons = set()
    for t in range(2 * k - 1):
        group = []
        for i in range(k):
            j = t - i
            if i <= j < k and j >= 0:
                a = [0] * k; a[i] += 1; a[j] += 1; group.append(tuple(a[1:-1]))
        groups.append(np.array(group, dtype=np.int64))
        for a, b in itertools.combinations(group, 2):
            z = primitive([x - y for x, y in zip(a, b)]); first = next(x for x in z if x)
            if first < 0:
                z = tuple(-x for x in z)
            comparisons.add(z)
    return d, groups, np.array(sorted(comparisons), dtype=np.int64)


def cross_directions(C, d, prefix=()):
    """All primitive directions orthogonal to a (d-1)-subset of comparison normals starting with `prefix`
    (exact: Python integers in object arrays)."""
    dirs = set(); n = len(C); need = d - 1 - len(prefix)
    if need < 0:
        return dirs
    rest = itertools.combinations(range(prefix[-1] + 1, n), need) if prefix else itertools.combinations(range(n), need)
    while True:
        batch = [tuple(prefix) + c for c in itertools.islice(rest, 20000)]
        if not batch:
            break
        A = C[np.array(batch)].astype(object)
        coeff = {(): np.ones(len(batch), dtype=object)}
        for p in range(1, d):
            nxt = {}
            for cols in itertools.combinations(range(d), p):
                z = np.zeros(len(batch), dtype=object)
                for m, j in enumerate(cols):
                    z = z + (-1) ** (p - 1 - m) * coeff[cols[:m] + cols[m + 1:]] * A[:, p - 1, j]
                nxt[cols] = z
            coeff = nxt
        V = np.stack([(-1) ** j * coeff[tuple(i for i in range(d) if i != j)] for j in range(d)], axis=1)
        for row in V:
            if all(x == 0 for x in row):
                continue
            v = primitive(row); s0 = next(x for x in v if x)
            if s0 < 0:
                v = tuple(-x for x in v)
            dirs.add(v)
    return dirs


def prefixes(n, d):
    """Split the (d-1)-subsets of range(n) by their first two elements (for parallel runs)."""
    if d - 1 >= 2:
        return [(a, b) for a in range(n) for b in range(a + 1, n) if n - b - 1 >= d - 3]
    return [(a,) for a in range(n)]


def is_fan_ray(k, u, groups, d):
    uu = np.array(u[1:-1], dtype=object); active = []
    for G in groups:
        z = G.astype(object) @ uu; m = min(z); Gm = G[[x == m for x in z]]
        active.extend((Gm[1:] - Gm[0]).tolist())
    return d == 1 or (len(active) >= d - 1 and int_rank(active, stop=d) == d - 1)


def rays_from_directions(k, dirs, groups, d):
    out = set()
    for v in dirs:
        for sgn in (1, -1):
            u = (0,) + tuple(sgn * x for x in v) + (0,)
            if is_fan_ray(k, u, groups, d):
                out.add(u)
    return out


def line_enumeration(k):
    """All rays of the fan of hat-P(k), k >= 3, by line enumeration. Returns (rays, number of candidate lines)."""
    d, groups, C = fan_setup(k)
    if d == 1:
        dirs = {(1,)}
    else:
        dirs = cross_directions(C, d)
    return rays_from_directions(k, dirs, groups, d), len(dirs)


def normalize_ray(u):
    """Representative with u(0) = u(k-1) = 0, primitive integer (matches fan_walk.normalize)."""
    k = len(u); w0 = u[0]; dd = u[-1] - w0
    return primitive([(k - 1) * (u[i] - w0) - dd * i for i in range(k)])


# ---------------- profiles and exact facet lists ----------------
def profiles(k):
    P = set()
    for s in itertools.permutations(range(k)):
        h = [0] * (2 * k - 1)
        for i, j in enumerate(s):
            h[i + j] += 1
        P.add(tuple(h))
    return sorted(P)


def canonical_inequality(k, b, g):
    """Canonical primitive form (b', g') of b + g.x >= 0 modulo the affine hull sum x = k, sum t x_t = k(k-1):
    g' has g'_{2k-3} = g'_{2k-2} = 0."""
    slope = F(g[-1]) - F(g[-2]); a = F(g[-1]) - slope * (2 * k - 2)
    return primitive([F(b) + a * k + slope * k * (k - 1)] + [F(v) - a - slope * t for t, v in enumerate(g)])


def ray_inequality(u):
    """The potential inequality of u as b + g.x >= 0: 2 sum u - sum_t v(t) x_t >= 0, in canonical form."""
    k = len(u)
    v = [min(u[i] + u[t - i] for i in range(k) if 0 <= t - i < k) for t in range(2 * k - 1)]
    return canonical_inequality(k, 2 * sum(u), [-x for x in v])


def coordinate_inequality(k, t):
    g = [0] * (2 * k - 1); g[t] = 1
    return canonical_inequality(k, 0, g)


def facet_audit(k, candidates, P=None):
    """Keep the candidate inequalities (b, g) that are valid on all profiles and tight on an affinely independent
    set of 2k - 3 profiles. Exact integer arithmetic."""
    P = P or profiles(k); V = np.array(P, dtype=object); d = 2 * k - 3; good = set()
    for h in sorted(set(candidates)):
        z = np.array(h[1:], dtype=object); slack = h[0] + V.dot(z)
        if min(slack) < 0:
            raise RuntimeError(f"inequality {h} is not valid")
        face = [P[i] for i in range(len(P)) if slack[i] == 0]
        if face and int_rank([[a - b for a, b in zip(f, face[0])] for f in face[1:]], stop=d - 1) == d - 1:
            good.add(h)
    return good


def exact_facets(k):
    """Exact facet list of hat-P(k) (canonical inequalities) from line enumeration plus the coordinates."""
    rays = set() if k <= 2 else line_enumeration(k)[0]
    H = {ray_inequality(u) for u in rays} | {coordinate_inequality(k, t) for t in range(2 * k - 1)}
    return facet_audit(k, H), rays


def profile_is_facet(u, sign=1, P=None):
    """Facet test by enumeration of all profiles (independent of the tight-graph test; small k only)."""
    k = len(u); uu = [sign * x for x in u]; P = P or profiles(k)
    v = [min(uu[i] + uu[t - i] for i in range(max(0, t - k + 1), min(t, k - 1) + 1)) for t in range(2 * k - 1)]
    rhs = 2 * sum(uu); lhs = [sum(a * b for a, b in zip(p, v)) for p in P]
    assert max(lhs) <= rhs
    T = [p for p, x in zip(P, lhs) if x == rhs]
    if len(T) == len(P) or not T:
        return False
    return int_rank([[a - b for a, b in zip(t, T[0])] for t in T[1:]], stop=2 * k - 4) == 2 * k - 4


# ---------------- floating-point double description (optional dependency pycddlib) ----------------
def cdd_facets_verified(k, P=None):
    """Facets of hat-P(k) by cddlib from all profiles (floating point); every row is rationalized to a primitive
    integer row (b, g) and verified exactly: valid on all profiles and tight on an affinely independent set of 2k-3
    profiles. Returns the set of canonical inequalities."""
    import cdd
    P = P or profiles(k)
    mat = cdd.matrix_from_array([[1] + list(r) for r in P], rep_type=cdd.RepType.GENERATOR)
    ineq = cdd.copy_inequalities(cdd.polyhedron_from_matrix(mat))
    rows = [r for i, r in enumerate(ineq.array) if i not in ineq.lin_set]
    V = np.array(P, dtype=object); out = set()
    for r in rows:
        q = [F(x).limit_denominator(10**6) for x in r]
        z = primitive(q)
        slack = z[0] + V.dot(np.array(z[1:], dtype=object))
        if min(slack) < 0:
            raise RuntimeError("cddlib row is not valid after rationalization")
        T = [P[i] for i in range(len(P)) if slack[i] == 0]
        if int_rank([[a - b for a, b in zip(t, T[0])] for t in T[1:]], stop=2 * k - 4) != 2 * k - 4:
            raise RuntimeError("cddlib row is not a facet after rationalization")
        out.add(canonical_inequality(k, z[0], z[1:]))
    if len(out) != len(rows):
        raise RuntimeError("duplicate rows")
    return out
