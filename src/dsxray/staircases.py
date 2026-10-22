"""Step potentials as staircases.

A composition g = (g_0, ..., g_r) of k gives the staircase u(g)(i) = l for i in block l (block l has length g_l);
every signed step potential is, up to sign and an affine term, the staircase of some composition with sign +.
The block ends are E_l = g_0 + ... + g_l - 1 and the block starts S_l = E_{l-1} + 1 (S_0 = 0).

This module provides the compositions, the involution g -> g* (composition of the complementary step set, mirrored),
the sufficient criterion S4, the classification rule for three steps, a certified sweep over all compositions of k,
and brute-force checks of every step in the proof of the three-step theorem.
"""
import itertools
from collections import Counter

from .facet_test import (allowed_edges, classify, perfect_matching, tight_edges)


def comps(k):
    """All compositions of k (as tuples), including the trivial ones."""
    for cuts in itertools.product([0, 1], repeat=k - 1):
        g = []; c = 1
        for b in cuts:
            if b:
                g.append(c); c = 1
            else:
                c += 1
        g.append(c)
        yield tuple(g)


def staircases(k):
    """Compositions of k with 1 <= r <= k-2 steps (r = number of parts - 1)."""
    return (g for g in comps(k) if 1 <= len(g) - 1 <= k - 2)


def staircases_r(k, r):
    """Compositions of k with exactly r steps (r + 1 parts), in lexicographic order of the cut points."""
    for cuts in itertools.combinations(range(1, k), r):
        b = (0,) + cuts + (k,)
        yield tuple(b[i + 1] - b[i] for i in range(r + 1))


def u_of(g):
    return [l for l, L in enumerate(g) for _ in range(L)]


def ends(g):
    E, s = [], -1
    for x in g:
        s += x; E.append(s)
    return E


def starts(g):
    E = ends(g)
    return [0] + [e + 1 for e in E[:-1]]


def maximizer_pattern(g):
    """For each level s = 1..2r-1, the block pairs (p, s-p) attaining M(s) = max E_p + E_{s-p}."""
    r = len(g) - 1; E = ends(g); pat = []
    for s in range(1, 2 * r):
        vals = {(p, s - p): E[p] + E[s - p] for p in range(max(0, s - r), s // 2 + 1)}
        m = max(vals.values())
        pat.append(tuple(pq for pq, v in sorted(vals.items()) if v == m))
    return tuple(pat)


def triple_rule(g):
    """Three-step theorem: u(g) with g = (g0, g1, g2, g3) is a facet unless g1, g2, g3 are pairwise distinct and g1 < g3."""
    _, a, b, c = g
    return not (len({a, b, c}) == 3 and a < c)


def J_of(g):
    s = 0; J = set()
    for x in g[:-1]:
        s += x; J.add(s)
    return J


def g_of(J, k):
    pts = [0] + sorted(J) + [k]
    return tuple(pts[i + 1] - pts[i] for i in range(len(pts) - 1))


def star(g):
    """g -> g*: the composition of {k - j : j in the complement of J(g)}."""
    k = sum(g); J = J_of(g)
    return g_of({k - j for j in range(1, k) if j not in J}, k)


def S4(g):
    """Sufficient criterion for a facet (conditions (U), (R) and g_0 >= 2)."""
    r = len(g) - 1
    if r < 1 or g[0] < 2:
        return False
    E = ends(g)
    M = lambda s: max(E[p] + E[s - p] for p in range(max(0, s - r), min(s, r) + 1))
    U = all(E[p] + E[s - p] <= E[0] + E[s] for s in range(r) for p in range(s + 1))
    R = all(E[l - 1] + E[l] >= M(2 * l - 1) for l in range(1, r + 1))
    return U and R


def total_support(g):
    u = u_of(g); k = len(u); E = tight_edges(u); m = perfect_matching(k, E)
    return m is not None and len(allowed_edges(k, E, m)) == len(E)


def classify_staircase(g):
    res = classify(u_of(g))
    res.update(g=list(g), k=sum(g), r=len(g) - 1)
    if len(g) == 4:
        res["triple_rule"] = triple_rule(g)
    return res


def summarize_classifications(rows):
    """Counts over certified classifications of staircases: facets, total support, facet == total support,
    certificates, rays, distinct maximizer patterns and patterns with both outcomes."""
    s = Counter(); certs = Counter(); pat = {}
    for r in rows:
        s["items"] += 1; s["facets"] += r["facet"]
        ts = bool(r.get("total_support")); s["total_support"] += ts
        s["facet_equals_total_support"] += (r["facet"] == ts); s["rays"] += r["ray"]
        if not r["facet"]:
            certs[r["cert"]] += 1
        if "triple_rule" in r:
            s["three_step"] += 1; s["three_step_rule_agrees"] += (r["facet"] == r["triple_rule"])
        key = (r["r"], maximizer_pattern(tuple(r["g"])))
        pat[key] = pat.get(key, 0) | (1 if r["facet"] else 2)
    out = dict(s)
    out["nonfacet_certificates"] = dict(sorted(certs.items()))
    out["maximizer_patterns"] = len(pat)
    out["patterns_with_both_outcomes"] = sum(1 for v in pat.values() if v == 3)
    return out


def nonincreasing_tails(k):
    """Staircases g = (g_0, g_1, ..., g_r) of k with g_1 >= g_2 >= ... >= g_r and 1 <= r <= k-2."""
    def dec_parts(n, r, maxpart):
        if r == 0:
            if n == 0:
                yield ()
            return
        for first in range(min(n - (r - 1), maxpart), 0, -1):
            for rest in dec_parts(n - first, r - 1, first):
                yield (first,) + rest
    for g0 in range(1, k):
        for r in range(1, min(k - 2, k - g0) + 1):
            for h in dec_parts(k - g0, r, k):
                yield (g0,) + h


# ---------------- three-step theorem: brute-force checks of every proof step ----------------
def _parts(g):
    E = [sum(g[:l + 1]) - 1 for l in range(4)]; S = [0] + [e + 1 for e in E[:-1]]
    L = [l for l, x in enumerate(g) for _ in range(x)]
    return E, S, L


def three_step_fails(g):
    return len({g[1], g[2], g[3]}) == 3 and g[1] < g[3]


def check_triangles(g, T):
    """Closed-form sizes m(p, q) of the corner triangles of the tight graph."""
    g0, g1, g2, g3 = g; E, S, L = _parts(g)
    pred = {(1, 2): min(g1, g2), (1, 3): min(g1, g3), (2, 3): min(g2, g3), (2, 2): min(g1 + g2 - g3, g2), (1, 1): g1, (3, 3): g3}
    tri = lambda m: m * (m + 1) // 2 if m > 0 else 0
    return all(sum(1 for (i, j) in T if L[i] == p and L[j] == q) == tri(m) for (p, q), m in pred.items())


def check_failing(g, T):
    """Failing cases: the stated neighbor sets of the starts S_2 and S_3."""
    g0, g1, g2, g3 = g; E, S, L = _parts(g); k = sum(g)
    nb = lambda i: {j for j in range(k) if (i, j) in T}
    if g2 < g1 < g3:
        return nb(S[3]) == {E[3]} and g3 >= 3
    if g1 < g2 < g3:
        return nb(S[2]) == {E[0], E[3]} and nb(S[3]) == {E[0], E[3]} and g3 >= 3
    return nb(S[2]) == {E[0]}                                    # g1 < g3 < g2


def involution(g):
    """The tight involution used in the facet cases (I), (II), (III)."""
    g0, g1, g2, g3 = g; E, S, L = _parts(g); idx = lambda l, a: E[l] - a; M = {}
    for a in range(g0):
        M[idx(0, a)] = idx(0, a)
    if g1 >= g3:
        for l in (1, 2, 3):
            for a in range(g[l]):
                M[idx(l, a)] = idx(l, g[l] - 1 - a)
    elif g1 == g2:
        for a in range(g1):
            M[idx(1, a)] = idx(2, g1 - 1 - a); M[idx(2, a)] = idx(1, g1 - 1 - a)
        for a in range(g3):
            M[idx(3, a)] = idx(3, g3 - 1 - a)
    else:
        for a in range(g1):
            M[idx(1, a)] = idx(1, g1 - 1 - a)
        for a in range(g2):
            M[idx(2, a)] = idx(3, g2 - 1 - a); M[idx(3, a)] = idx(2, g2 - 1 - a)
    return M


def check_facet_steps(g, T):
    """Facet cases: the involution is a tight permutation, its switching digraph is strongly connected, and both
    step classes (within-block and boundary) are connected by square relations."""
    k = sum(g); E, S, L = _parts(g); M = involution(g)
    if sorted(M.values()) != list(range(k)) or not all((i, M[i]) in T for i in range(k)):
        return False
    inv = {v: i for i, v in M.items()}
    succ = {i: {inv[j] for j in range(k) if (i, j) in T} for i in range(k)}

    def reach(s):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in succ[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        return seen
    if not (reach(E[0]) == set(range(k)) and all(E[0] in reach(i) for i in range(k))):
        return False
    par = list(range(k - 1))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(k - 1):
        for j in range(k - 1):
            if (i, j + 1) in T and (i + 1, j) in T:
                par[f(i)] = f(j)
    return (len({f(i) for i in range(k - 1) if L[i + 1] == L[i]}) <= 1
            and len({f(i) for i in range(k - 1) if L[i + 1] != L[i]}) == 1)


def three_step_proof_checks(kmin, kmax, theorem_check=True):
    """Counts (checked, passed) for each proof step and for the statement, over all three-step staircases."""
    out = {"triangles": [0, 0], "failing": [0, 0], "facet_steps": [0, 0], "theorem": [0, 0]}
    for k in range(max(5, kmin), kmax + 1):
        for x, y, z in itertools.combinations(range(1, k), 3):
            g = (x, y - x, z - y, k - z); T = tight_edges(u_of(g))
            out["triangles"][0] += 1; out["triangles"][1] += check_triangles(g, T)
            fl = three_step_fails(g)
            key = "failing" if fl else "facet_steps"
            out[key][0] += 1; out[key][1] += (check_failing(g, T) if fl else check_facet_steps(g, T))
            if theorem_check:
                out["theorem"][0] += 1; out["theorem"][1] += (classify(u_of(g))["facet"] == (not fl))
    return out
