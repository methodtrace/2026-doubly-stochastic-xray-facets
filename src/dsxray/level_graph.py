"""Level structure of a staircase: block ends, the sequence M(s), attaining pairs, the classes of boundary steps,
statement (H+) (not used in the paper), and the rigid-corner certificate of conj:E.

For a composition g with block ends E_0 < ... < E_r put M(s) = max_{p+q=s} E_p + E_q (0 <= s <= 2r). A pair (p, q)
attains level s = p + q if E_p + E_q = M(s). By the corner lemma, (i, j) with i in block p and j in block q is tight
iff i + j > M(p + q - 1) (always when p = q = 0), and the tight pairs on the anti-diagonal M(s) are exactly the corners
(E_c, E_d) of the attaining pairs of level s.
"""
import itertools
from collections import Counter

from .facet_test import allowed_edges, perfect_matching, rank_mod, tight_edges, v_of, z_rows
from .staircases import ends, staircases, u_of


def Mfun(E):
    r = len(E) - 1
    return [max(E[p] + E[s - p] for p in range(max(0, s - r), min(s, r) + 1)) for s in range(2 * r + 1)]


def attainers(E, s):
    """Unordered attaining pairs (c, d), c <= d, of level s."""
    r = len(E) - 1
    L = [(c, s - c) for c in range(max(0, s - r), min(r, s) + 1) if c <= s - c]
    m = max(E[c] + E[d] for c, d in L)
    return [(c, d) for c, d in L if E[c] + E[d] == m]


def attain_matrix(E, M=None):
    r = len(E) - 1; M = M or Mfun(E)
    return [[E[p] + E[q] == M[p + q] for q in range(r + 1)] for p in range(r + 1)]


def boundary_classes(E, att=None):
    """Classes of the boundary steps {0..r-1}: p ~ p' when both occur in attaining pairs of one common level."""
    r = len(E) - 1; att = att or attain_matrix(E)
    par = list(range(r))

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for s in range(2 * r + 1):
        mem = [p for p in range(max(0, s - r), min(s, r) + 1) if att[p][s - p] and p < r]
        for x in mem[1:]:
            par[f(x)] = f(mem[0])
    cl = {}
    for p in range(r):
        cl.setdefault(f(p), []).append(p)
    return sorted(cl.values())


def gstar_tight_graph_reflects(g):
    """The tight graph of g* is the reflection (i, j) -> (k-1-i, k-1-j) of the tight graph of g."""
    from .staircases import star
    k = sum(g)
    T = tight_edges(u_of(g)); Ts = tight_edges(u_of(star(g)))
    return Ts == {(k - 1 - i, k - 1 - j) for i, j in T}


def start_neighbors_are_attaining_ends(g):
    """Lemma on block starts: the tight neighbors of S_{q+1} = E_q + 1 are {E_l : (q, l) attains M(q + l)}."""
    k = sum(g); r = len(g) - 1; E = ends(g); M = Mfun(E); T = tight_edges(u_of(g))
    for q in range(r):
        nb = sorted(j for j in range(k) if (E[q] + 1, j) in T)
        if nb != sorted(E[l] for l in range(r + 1) if E[q] + E[l] == M[q + l]):
            return False
    return True


def last_antidiagonal_corners(g):
    """Tight pairs on the anti-diagonal M(s) are exactly the corners of the attaining block pairs of level s."""
    k = sum(g); r = len(g) - 1; u = u_of(g); v = v_of(u); E = ends(g)
    for s in range(2 * r + 1):
        at = attainers(E, s); t = E[at[0][0]] + E[at[0][1]]
        tight = {(i, t - i) for i in range(max(0, t - k + 1), min(k, t + 1)) if u[i] + u[t - i] == v[t]}
        corners = {(E[c], E[d]) for c, d in at} | {(E[d], E[c]) for c, d in at}
        if tight != corners:
            return False
    return True


# ---------------- square-relation graphs ----------------
def _ncomponents(nodes, edges):
    par = {x: x for x in nodes}

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in edges:
        par[f(a)] = f(b)
    return len({f(x) for x in nodes})


def square_relation_analysis(g):
    """Total support, connectivity of the boundary-step graph under the square relation, connectivity of the
    within-block steps, the ray property, and existence of a tight permutation."""
    k = sum(g); r = len(g) - 1; E = ends(g); M = Mfun(E)
    bnodes = list(range(r)); bedges = [(p, q) for p in range(r) for q in range(p + 1, r) if E[p] + E[q] == M[p + q]]
    mtri = lambda p, q: None if p + q == 0 else max(0, E[p] + E[q] - M[p + q - 1])
    wnodes = [(p, a) for p in range(r + 1) for a in range(1, g[p])]
    wedges = [((p, a), (q, b)) for (p, a), (q, b) in itertools.combinations(wnodes, 2)
              if (p == 0 and q == 0) or a + b <= mtri(p, q)]
    u = u_of(g); Et = tight_edges(u); mt = perfect_matching(k, Et)
    ts = mt is not None and len(allowed_edges(k, Et, mt)) == len(Et)
    zd = k - rank_mod(z_rows(k, Et))
    return {"total_support": ts, "boundary_connected": _ncomponents(bnodes, bedges) == 1,
            "within_connected": (_ncomponents(wnodes, wedges) == 1) if wnodes else True,
            "ray": zd == 3, "tight_permutation": mt is not None}


# ---------------- statement (H+) ----------------
def hplus_verdict(g):
    """Statement (H+) for a composition: 'connected' (one class), 'deficient' (some class K has |A(K)| < |K|),
    'tight_set' (a balanced class with a tight edge from an end of A(K) to an index outside the starts of K), or
    'FAIL' (several classes and neither alternative). Tightness is decided by the corner lemma."""
    r = len(g) - 1
    E = [e - 1 for e in itertools.accumulate(g)]; S = [E[p] - g[p] + 1 for p in range(r + 1)]
    M = Mfun(E); att = attain_matrix(E, M)
    comps = boundary_classes(E, att)
    if len(comps) == 1:
        return "connected", None

    def tight_with_outside(l, X):
        for p in range(r + 1):
            lo = S[p] if (p, l) != (0, 0) else 0
            if p + l > 0:
                lo = max(lo, M[p + l - 1] + 1 - E[l])
            for j in range(lo, E[p] + 1):
                if j not in X:
                    return True
        return False
    for K in comps:
        A = {l for q in K for l in range(r + 1) if att[q][l]}
        if len(A) < len(K):
            return "deficient", K
        if len(A) == len(K):
            X = {S[q + 1] for q in K}
            if any(tight_with_outside(l, X) for l in A):
                return "tight_set", K
    return "FAIL", comps


def hplus_counts(k):
    C = Counter(); fails = []
    for g in staircases(k):
        res, K = hplus_verdict(g); C[res] += 1
        if res == "FAIL" and len(fails) < 20:
            fails.append(list(g))
    return {"k": k, "counts": dict(sorted(C.items())), "fails": fails}


# ---------------- rigid-corner certificate ----------------
def rigid_corner(g):
    """'total_support', 'no_pm' (no tight permutation), 'certified' (some level s has a unique attaining pair (c, d)
    whose corner (E_c, E_d) lies in no tight permutation), or 'FAIL'."""
    k = sum(g); r = len(g) - 1; T = tight_edges(u_of(g))
    m = perfect_matching(k, T)
    if m is None:
        return "no_pm"
    A = allowed_edges(k, T, m)
    if len(A) == len(T):
        return "total_support"
    E = ends(g)
    for b in range(2 * r + 1):
        at = attainers(E, b)
        if len(at) == 1:
            c, d = at[0]
            if (E[c], E[d]) not in A:
                return "certified"
    return "FAIL"


def first_single_start_partner(g):
    """None if no block start S_l (l >= 1) has a unique tight neighbor E_q with 0 < q < r (taking the first such
    start); otherwise whether some p with E_p != S_l makes (p, q) the unique attaining pair of level p + q."""
    k = sum(g); r = len(g) - 1; u = u_of(g); v = v_of(u); E = ends(g); S = [0] + [e + 1 for e in E[:-1]]
    for l in range(1, r + 1):
        s = S[l]; N = [j for j in range(k) if u[s] + u[j] == v[s + j]]
        if len(N) != 1:
            continue
        q = u[N[0]]
        if q in (0, r):
            return None
        return any(E[p] != s and len(attainers(E, p + q)) == 1 and tuple(sorted((p, q))) == attainers(E, p + q)[0]
                   for p in range(r + 1))
    return None


def several_classes_without_deficient_class(g):
    """(the boundary steps form several classes, no class K has |A(K)| < |K|)."""
    r = len(g) - 1; E = ends(g); att = attain_matrix(E)
    comps = boundary_classes(E, att)
    if len(comps) == 1:
        return False, False
    deficient = any(len({l for q in K for l in range(r + 1) if att[q][l]}) < len(K) for K in comps)
    return True, not deficient
