"""Certified facet test for the potential inequalities of the X-ray polytope hat-P(k).

A potential u on {0..k-1} gives the min-plus inequality  sum_t v(t) x_t <= 2 sum_i u(i),  v(t) = min_{i+j=t} u(i)+u(j),
valid on hat-P(k) (laws of i + sigma(i), sigma a permutation, in count coordinates). The tight graph E(u) is the set of
pairs (i, j) with u(i) + u(j) = v(i+j); a permutation attains the bound iff it lies in E(u). The allowed graph is the
set of tight edges that lie in some tight permutation (Dulmage-Mendelsohn: same strong component of the switching
digraph). With M the row/column incidence matrix of the allowed edges and A the anti-diagonal (X-ray) rows,

    dim(face) = rank [M; A] - rank M,       facet iff dim = 2k - 4 and the inequality is not an equation.

Every outcome of `classify` carries a certificate:
  "empty"    no tight permutation (Hall failure): the face is empty;
  "coords"   two anti-diagonals carry no allowed edge: the face lies in two independent coordinate hyperplanes;
  "equation" every edge is tight: the inequality holds with equality on hat-P(k);
  "rank"     face dimension 2k - 4 from a rank modulo the prime 2^31 - 1. M is totally unimodular, so its rank is
             the same modulo every prime, and a modular rank never exceeds the rational rank; so reaching 2k - 4
             certifies a facet over Q;
  "exact"    the modular rank fell short and the dimension was recomputed in exact rational arithmetic.
The ray test (dim Z(u) = 3) is certified the same way.

`reference_is_facet` is a slower, independent implementation (allowed edges by deleting each edge and testing for a
perfect matching, ranks modulo 2^61 - 1 and in rational arithmetic); it shares no code with `classify`.
"""
import sys
from fractions import Fraction

import numpy as np

P31 = 2147483647
P61 = (1 << 61) - 1


# ---------------- tight graph ----------------
def v_of(u):
    k = len(u)
    return [min(u[i] + u[t - i] for i in range(max(0, t - k + 1), min(t, k - 1) + 1)) for t in range(2 * k - 1)]


def tight_edges(u):
    k = len(u); v = v_of(u)
    return {(i, j) for i in range(k) for j in range(k) if u[i] + u[j] == v[i + j]}


def perfect_matching(k, E):
    """A perfect matching of the bipartite graph E on rows x columns {0..k-1}, as a dict column -> row, or None."""
    adj = [[j for j in range(k) if (i, j) in E] for i in range(k)]
    match = {}

    def aug(i, seen):
        for j in adj[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in match or aug(match[j], seen):
                match[j] = i
                return True
        return False
    sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
    for i in range(k):
        if not aug(i, set()):
            return None
    return match


def allowed_edges(k, E, match):
    """Edges of E that lie in some perfect matching, given one perfect matching (column -> row)."""
    rowof = match
    adj = {i: [j for j in range(k) if (i, j) in E] for i in range(k)}
    succ = {i: {rowof[j] for j in adj[i]} for i in range(k)}
    reach = []
    for s in range(k):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in succ[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        reach.append(seen)
    return {(i, j) for (i, j) in E if i == rowof[j] or (rowof[j] in reach[i] and i in reach[rowof[j]])}


def has_total_support(u):
    k = len(u); E = tight_edges(u); m = perfect_matching(k, E)
    return m is not None and len(allowed_edges(k, E, m)) == len(E)


# ---------------- ranks ----------------
def rank_mod(rows, p=P31):
    """Rank modulo the prime p (numpy, int64; p < 2^31 keeps every product below 2^62)."""
    if not rows:
        return 0
    A = np.array(rows, dtype=np.int64) % p
    if A.shape[1] > A.shape[0]:
        A = A.T.copy()
    r = 0; R, C = A.shape
    for c in range(C):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]; A[[r, piv]] = A[[piv, r]]
        A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        f = A[:, c].copy(); f[r] = 0
        A = (A - (f[:, None] * A[r][None, :]) % p) % p
        r += 1
        if r == R:
            break
    return r


def rank_exact(rows):
    """Rank over Q (fractions.Fraction)."""
    rows = [[Fraction(x) for x in r] for r in rows]
    rank = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((r for r in range(rank, len(rows)) if rows[r][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(len(rows)):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c] / rows[rank][c]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rank


def face_rows(k, edges):
    Ep = sorted(edges)
    M = [[1 if e[0] == i else 0 for e in Ep] for i in range(k)] + [[1 if e[1] == j else 0 for e in Ep] for j in range(k)]
    A = [[1 if e[0] + e[1] == t else 0 for e in Ep] for t in range(2 * k - 1)]
    return M, A


def z_rows(k, E):
    """Rows of the system defining Z(u) = {z : z(i) + z(j) is constant on the tight pairs of each anti-diagonal}."""
    rows = []
    for t in range(2 * k - 1):
        Pt = [(i, t - i) for i in range(k) if 0 <= t - i < k and (i, t - i) in E]
        for (i, j) in Pt[1:]:
            r = [0] * k; r[i] += 1; r[j] += 1; r[Pt[0][0]] -= 1; r[Pt[0][1]] -= 1
            rows.append(r)
    return rows


def z_dim(u):
    """dim Z(u) in exact arithmetic; for a non-affine u it equals 3 exactly when u spans a ray of the fan."""
    k = len(u); rows = z_rows(k, tight_edges(u))
    return k - (rank_exact(rows) if rows else 0)


def is_ray(u):
    return z_dim(u) == 3


# ---------------- certified classification ----------------
def classify(u):
    """Certified classification of the inequality of a (non-affine) potential u. Returns a dict with keys
    ray, zdim, facet, cert, total_support, and (when a tight permutation exists) lost_edges, empty_antidiagonals."""
    k = len(u); E = tight_edges(u); out = {}
    zr = z_rows(k, E); zd = k - rank_mod(zr)
    if zd != 3:
        zd = k - rank_exact(zr)
    out["ray"] = (zd == 3); out["zdim"] = zd
    m = perfect_matching(k, E)
    if m is None:
        out.update(facet=False, cert="empty", total_support=False, empty_antidiagonals=None)
        return out
    Ep = allowed_edges(k, E, m)
    out["total_support"] = (len(Ep) == len(E)); out["lost_edges"] = len(E) - len(Ep)
    empty = [t for t in range(2 * k - 1) if not any((i, t - i) in Ep for i in range(k) if 0 <= t - i < k)]
    out["empty_antidiagonals"] = len(empty)
    if len(empty) >= 2:
        out.update(facet=False, cert="coords")
        return out
    if len(Ep) == k * k:
        out.update(facet=False, cert="equation")
        return out
    M, A = face_rows(k, Ep)
    d = rank_mod(M + A) - rank_mod(M)
    if d == 2 * k - 4:
        out.update(facet=True, cert="rank")
        return out
    d = rank_exact(M + A) - rank_exact(M)
    out.update(facet=(d == 2 * k - 4), cert="exact", face_dim=d)
    return out


def classify_coordinate(k, t):
    """The coordinate inequality x_t >= 0: its face is spanned by the permutations avoiding anti-diagonal t."""
    E = {(i, j) for i in range(k) for j in range(k) if i + j != t}
    m = perfect_matching(k, E)
    if m is None:
        return {"t": t, "facet": False, "cert": "empty"}
    Ep = allowed_edges(k, E, m); M, A = face_rows(k, Ep)
    d = rank_mod(M + A) - rank_mod(M)
    if d != 2 * k - 4:
        d = rank_exact(M + A) - rank_exact(M)
    return {"t": t, "facet": d == 2 * k - 4, "cert": "rank" if d == 2 * k - 4 else "exact", "face_dim": d}


def is_facet(u, sign=1):
    return classify([sign * x for x in u])["facet"]


# ---------------- independent reference implementation ----------------
def _has_pm(k, E, skip_row=None, skip_col=None):
    adj = {i: [j for j in range(k) if (i, j) in E and j != skip_col] for i in range(k) if i != skip_row}
    match = {}

    def aug(i, seen):
        for j in adj[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in match or aug(match[j], seen):
                match[j] = i
                return True
        return False
    return all(aug(i, set()) for i in adj)


def matchable_edges(k, E):
    """Edges in some perfect matching, by deleting the edge's row and column and testing for a perfect matching."""
    if not _has_pm(k, E):
        return set()
    return {(i, j) for (i, j) in E if _has_pm(k, E, i, j)}


def _rank_mod_list(rows, p=P61):
    rows = [[x % p for x in r] for r in rows]
    rank, ncol = 0, len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((r for r in range(rank, len(rows)) if rows[r][c]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rows[rank] = [x * inv % p for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][c]:
                f = rows[r][c]; rows[r] = [(a - f * b) % p for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rank


def reference_face_dim(u, exact=False):
    """(dimension of the tight face, is_equation); dimension -1 if no permutation is tight."""
    k = len(u); Ep = sorted(matchable_edges(k, tight_edges(u)))
    if not Ep:
        return -1, False
    M, A = face_rows(k, Ep)
    rk = rank_exact if exact else _rank_mod_list
    return rk(M + A) - rk(M), len(Ep) == k * k


def reference_is_facet(u, sign=1):
    uu = [sign * x for x in u]; k = len(u)
    d, eq = reference_face_dim(uu)
    if eq:
        return False
    if d == 2 * k - 4:
        return True
    return reference_face_dim(uu, exact=True)[0] == 2 * k - 4
