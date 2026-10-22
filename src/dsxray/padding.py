"""End padding of facet rays and the decomposition of facets into cores.

A ray u (nondecreasing after an affine change) is stored as its pattern p: increments d(i) = u(i+1) - u(i), shifted
to minimum 0 and divided by their gcd; the height H is max p. padL(p) = (0, p) duplicates the left end value;
padR(p) = (p, H) is its mirror image (mir(p) = H - reversed p). With E(u) the tight graph:
  P_L(u): some m has d(m) = 0 and (0, m) in E(u);  P_R(u) = P_L(mirror of u).

Padding theorem: if u is a ray with connected tight graph and total support and P_L(u) holds, then padL(u) is a ray
with the same properties and gives a facet of hat-P(k+1); padL preserves P_L and P_R (and mirror statements). A core
is a facet ray that is neither padL nor padR of a facet ray one size smaller; every facet with k <= 12 is
padL^a padR^b of a unique core.
"""
from functools import reduce
from math import gcd

from .facet_test import allowed_edges, classify, perfect_matching, tight_edges


def increments(u):
    return [u[i + 1] - u[i] for i in range(len(u) - 1)]


def pattern(u):
    d = increments(u); m = min(d); d = [x - m for x in d]; g = reduce(gcd, d) or 1
    return tuple(x // g for x in d)


def from_pattern(p):
    u = [0]
    for x in p:
        u.append(u[-1] + x)
    return u


def is_step(p):
    return max(p) <= 1


def pad_left(p):
    return (0,) + tuple(p)


def pad_right(p):
    return tuple(p) + (max(p),)


def mirror(p):
    H = max(p)
    return tuple(H - x for x in reversed(p))


def normalized(q):
    return len(set(q)) > 1 and pattern(from_pattern(q)) == tuple(q)


def P_L(p):
    E = tight_edges(from_pattern(p))
    return any(p[m] == 0 and (0, m) in E for m in range(len(p)))


def P_R(p):
    H = max(p); k = len(p) + 1; E = tight_edges(from_pattern(p))
    return any(p[m] == H and (m + 1, k - 1) in E for m in range(len(p)))


def connected_with_total_support(u):
    """(total support, bipartite tight graph connected)."""
    k = len(u); E = tight_edges(u); m = perfect_matching(k, E)
    if m is None or len(allowed_edges(k, E, m)) != len(E):
        return False, None
    adj = {("r", i): set() for i in range(k)}; adj.update({("c", j): set() for j in range(k)})
    for i, j in E:
        adj[("r", i)].add(("c", j)); adj[("c", j)].add(("r", i))
    seen = {("r", 0)}; st = [("r", 0)]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    return True, len(seen) == 2 * k


def is_facet_pattern(p, sign=1):
    return classify([sign * x for x in from_pattern(p)])["facet"]


def strip_decompositions(p, facets_by_k):
    """All maximal ways of stripping p by inverse paddings while staying among facets: set of (a, b, core)."""
    def canL(q):
        return q[0] == 0 and normalized(q[1:]) and tuple(q[1:]) in facets_by_k.get(len(q), set())

    def canR(q):
        return q[-1] == max(q) and normalized(q[:-1]) and tuple(q[:-1]) in facets_by_k.get(len(q), set())
    seen = {}; st = [(tuple(p), 0, 0)]; res = set()
    while st:
        q, a, b = st.pop()
        if (a, b) in seen:
            continue
        seen[(a, b)] = q; moved = False
        if len(q) >= 2 and canL(q):
            st.append((q[1:], a + 1, b)); moved = True
        if len(q) >= 2 and canR(q):
            st.append((q[:-1], a, b + 1)); moved = True
        if not moved:
            res.add((a, b, q))
    return res


def core_table(facets_by_k, kmin, kmax):
    """facets_by_k: k -> set of facet patterns (length k-1). Returns per-k counts and the cores."""
    table = {}; cores = {}
    for k in range(kmin, kmax + 1):
        Fk = facets_by_k[k]; nonuniq = 0; cs = []
        for p in Fk:
            r = strip_decompositions(p, facets_by_k)
            if len({c for _, _, c in r}) > 1:
                nonuniq += 1
            if (0, 0, p) in r:
                cs.append(p)
        ns = [c for c in cs if not is_step(c)]
        nsF = [p for p in Fk if not is_step(p)]
        from_padding = sum(1 for p in nsF if p not in set(cs))
        table[k] = {"ray_facets": len(Fk), "cores": len(cs), "nonstep_facets": len(nsF), "nonstep_cores": len(ns),
                    "nonstep_from_padding": from_padding, "nonunique_core": nonuniq}
        cores[k] = sorted(cs)
    return table, cores


def small_facet_patterns(k):
    """Facet patterns for k = 3, 4 by the certified test over all patterns in {0,1,2}^(k-1)."""
    import itertools
    S = set()
    for p in itertools.product(range(3), repeat=k - 1):
        if min(p) != 0 or max(p) == 0 or pattern(from_pattern(p)) != p:
            continue
        if is_facet_pattern(p):
            S.add(p)
    return S


def padding_structure_ok(p):
    """Steps (i) and (ii) of the proof of the padding theorem for u = from_pattern(p), u' = padL(u), checked on the
    computed min-plus self-convolutions and tight graphs:
      v'(0) = v'(1) = 2u(0) and v'(t) = v(t-2) for 2 <= t <= 2k;
      E(u') = shift of E(u) plus (0,0), (0,1), (1,0) plus (0, m+2), (m+2, 0) for m in M,
      M = {m : d(m) = 0 and (0, m) in E(u)}."""
    from .facet_test import v_of
    u = from_pattern(p); k = len(u); up = from_pattern(pad_left(p))
    v = v_of(u); vp = v_of(up)
    ok = vp[0] == vp[1] == 2 * u[0] and all(vp[t] == v[t - 2] for t in range(2, 2 * k + 1))
    E = tight_edges(u); d = increments(u)
    M = [m for m in range(k - 1) if d[m] == 0 and (0, m) in E]
    pred = {(i + 1, j + 1) for i, j in E} | {(0, 0), (0, 1), (1, 0)} | {(0, m + 2) for m in M} | {(m + 2, 0) for m in M}
    return ok and tight_edges(up) == pred
