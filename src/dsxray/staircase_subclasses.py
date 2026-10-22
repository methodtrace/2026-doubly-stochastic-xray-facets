"""Infinite classes of staircases on which E-flat holds (facet implies total support), checked exactly.

Notation of sec:census (cell_enumeration.py): a tail t = (g_1..g_r), ends E_0 = 0, E_l = g_1 + ... + g_l (this module
uses E_0 = 0, i.e. g_0 = 1; the classes below do not depend on g_0), M(s) = max_{p+q=s} E_p + E_q, the attainment
graph A(q) = {l : E_q + E_l = M(q + l)} for rows q = 0..r-1, balanced sets K (|A(K)| = |K|), K* the largest one.

  unimodal tails (g_1 <= ... <= g_m >= ... >= g_r): every nonempty balanced K has 0 in A(K); the key step is that a
      row q with 0 not in A(q) has {q, q + 1} in A(q) (a concave majorant from a maximizing prefix average);
  rigid neighbour of the first balanced row: if c in A(min K) and level min K + c has a unique maximizing pair, it
      certifies E-flat; hence E-flat holds whenever every level has a unique maximizing pair, in particular when
      there are no additive ties E_a + E_b = E_c + E_d, e.g. for any ordering of superincreasing parts;
  superadditive prefixes (H_p + H_q <= H_{p+q} for p + q <= r - 1, H_j = g_1 + ... + g_j): every row sees column 0;
  bipartite component: every staircase with a tight permutation but without total support has a bipartite allowed
      component (the vector 1_X - 1_N(X), X the ends of the largest balanced set, lies in Y*(u) and is nonzero).
`selector_counterexample` is the tail (4,5,4,1,5,5,5,5,1,1) on which the smallest block of A* has no rigid partner.
"""
from fractions import Fraction

from . import exact_graphs as GR


def ends_of(tail):
    E = [0]
    for x in tail:
        E.append(E[-1] + x)
    return E


def attainment(tail):
    """(E, A, pairs, M): A[q] the attaining neighbours of row q, pairs[s] the maximizing unordered pairs at level s."""
    r = len(tail); E = ends_of(tail)
    A = [set() for _ in range(r)]; pairs = []; M = []
    for s in range(2 * r + 1):
        ch = [(E[a] + E[s - a], a, s - a) for a in range(r + 1) if a <= s - a <= r]
        m = max(x[0] for x in ch); p = [(a, b) for v, a, b in ch if v == m]
        M.append(m); pairs.append(p)
        for a, b in p:
            if a < r:
                A[a].add(b)
            if b < r:
                A[b].add(a)
    return E, A, pairs, M


def balanced_by_subsets(A):
    """(Hall holds, largest balanced set) by enumerating all row subsets."""
    r = len(A); largest = set()
    for mask in range(1, 1 << r):
        K = {q for q in range(r) if mask >> q & 1}
        AK = set().union(*(A[q] for q in K))
        if len(AK) < len(K):
            return False, None
        if len(AK) == len(K):
            largest |= K
    return True, largest


def balanced_by_matching(A):
    """Largest balanced set by a maximum matching (rows without an alternating path to a free column); None if Hall
    fails."""
    r = len(A); mc = {}

    def aug(q, seen):
        for c in sorted(A[q]):
            if c in seen:
                continue
            seen.add(c)
            if c not in mc or aug(mc[c], seen):
                mc[c] = q
                return True
        return False
    if not all(aug(q, set()) for q in range(r)):
        return None
    mr = {q: c for c, q in mc.items()}
    rc = set(range(r + 1)) - set(mc); rr = set(); todo = list(rc)
    while todo:
        c = todo.pop()
        for q in range(r):
            if c in A[q] and mr[q] != c and q not in rr:
                rr.add(q)
                if mr[q] not in rc:
                    rc.add(mr[q]); todo.append(mr[q])
    return set(range(r)) - rr


def unimodal(tail):
    dec = False
    for a, b in zip(tail, tail[1:]):
        if b < a:
            dec = True
        elif b > a and dec:
            return False
    return True


def check_unimodal(tail):
    """Checks the steps of the unimodal-tail argument; returns the verdict ('Hall_fails', 'strict_Hall' or
    'endpoint_certificate'). Raises AssertionError if a step fails."""
    r = len(tail); E, A, pairs, M = attainment(tail)
    avg = [Fraction(E[j], j) for j in range(1, r + 1)]
    mx = max(avg); maxi = [j for j in range(1, r + 1) if avg[j - 1] == mx]; h = min(maxi)
    Fm = [mx * j if j < h else Fraction(E[j]) for j in range(r + 1)]
    sl = [Fm[j + 1] - Fm[j] for j in range(r)]
    assert all(Fm[j] >= E[j] for j in range(r + 1)) and all(x >= y for x, y in zip(sl, sl[1:])), "majorant"
    for q in range(r):
        if 0 not in A[q]:
            assert max(maxi) <= q and {q, q + 1} <= A[q], "row implication"
    hall, K = balanced_by_subsets(A)
    if hall:
        for mask in range(1, 1 << r):
            KK = {q for q in range(r) if mask >> q & 1}
            AK = set().union(*(A[q] for q in KK))
            if len(AK) == len(KK):
                assert 0 in AK, "balanced set without column 0"
        if K:
            assert pairs[0] == [(0, 0)]
            return "endpoint_certificate"
        return "strict_Hall"
    return "Hall_fails"


def no_additive_ties(tail):
    r = len(tail); E = ends_of(tail)
    return all(len(v) == len(set(v)) for v in ([E[a] + E[s - a] for a in range(r + 1) if a <= s - a <= r]
                                                for s in range(2 * r + 1)))


def rigid_first_row(tail, check_subsets=False):
    """Verdict of the rigid-neighbour certificate: 'hall_fails', 'balanced_empty', 'rigid_certificate' or
    'first_row_no_rigid_neighbor'."""
    r = len(tail); E, A, pairs, M = attainment(tail)
    K = balanced_by_matching(A)
    if check_subsets:
        hall, K2 = balanced_by_subsets(A)
        assert (K is None) == (not hall) and (K is None or K == K2)
    if K is None:
        return "hall_fails"
    if not K:
        return "balanced_empty"
    AK = set().union(*(A[q] for q in K))
    B = {d for d in range(1, r + 1) if tail[d - 1] == 1 and d - 1 in K}
    q = min(K)
    ch = [c for c in sorted(A[q]) if len(pairs[q + c]) == 1]
    if not ch:
        return "first_row_no_rigid_neighbor"
    c = ch[0]
    assert c in AK and q not in B and pairs[q + c] == [(min(q, c), max(q, c))]
    return "rigid_certificate"


def superadditive_prefix(g):
    """g = (g_0, g_1, ..., g_r): H_j = g_1 + ... + g_j for j <= r - 1; H_p + H_q <= H_{p+q} whenever p + q <= r - 1."""
    pref = [0]
    for x in g[1:-1]:
        pref.append(pref[-1] + x)
    return all(pref[p] + pref[q] <= pref[p + q] for p in range(len(pref)) for q in range(len(pref) - p))


def tail_structure(g):
    """For a composition g = (g_0..g_r): staircase ends, attainment neighbours, attaining pairs and K* (None if
    Hall fails), with E_l the last index of block l."""
    r = len(g) - 1; ends = []; acc = -1
    for x in g:
        acc += x; ends.append(acc)
    grp = [[(p, s - p) for p in range(max(0, s - r), s // 2 + 1)] for s in range(2 * r + 1)]
    mx = [max(ends[p] + ends[q] for p, q in G) for G in grp]
    att = [[(p, q) for p, q in G if ends[p] + ends[q] == mx[s]] for s, G in enumerate(grp)]
    nbs = [{l for l in range(r + 1) if ends[q] + ends[l] == mx[q + l]} for q in range(r)]
    hall, K = balanced_by_subsets(nbs)
    return ends, nbs, att, (sorted(K) if hall else None)


def staircase_u(g):
    return [lev for lev, n in enumerate(g) for _ in range(n)]


def bipartite_component_check(g):
    """For a staircase: None if no tight permutation, 'total_support', or (after checking that X = ends of K* + 1,
    N = N(X) give a nonzero y = 1_X - 1_N in Y*(u) and that a bipartite allowed component exists)
    'bipartite_component'; also checks the superadditive-prefix statements."""
    u = staircase_u(g); k = len(u)
    _, tight = GR.tight_graph(u)
    allowed = GR.allowed_edges(k, tight)
    ends, nbs, att, K = tail_structure(g)
    sup = superadditive_prefix(g)
    if sup:
        assert all(0 in nb for nb in nbs)
        if K:
            assert any(0 in nbs[q] for q in K) and att[0] == [(0, 0)]
    if allowed is None:
        return None
    if allowed == tight:
        return "total_support"
    assert K, "K* empty without total support"
    X = {ends[q] + 1 for q in K}; N = {ends[l] for q in K for l in nbs[q]}
    assert len(X) == len(N) and X != N
    y = [int(i in X) - int(i in N) for i in range(k)]
    assert any(y) and all(y[i] + y[j] == 0 for i, j in allowed)
    assert any(c["bipartite"] for c in GR.components(k, allowed))
    if sup:
        assert not any(i + j == 2 * ends[0] for i, j in allowed)
    return "bipartite_component"


EXAMPLES = {
    "height_two_k17": ([0, 0, 2, 2, 0, 0, 2, 2, 1, 1, 1, 0, 0, 2, 2, 0, 0], (6, 1, 28, 4)),
    "total_support_bipartite_k17": ([0, 0, 4, 4, 0, 0, 4, 4, 1, 2, 3, 0, 0, 4, 4, 0, 0], (6, 1, 28, 0)),
    "facet_without_total_support_k19": ([int(i in {2, 4, 5, 7, 9, 10, 11, 12, 13}) for i in range(19)], (4, 1, 34, 19)),
}


def example_diagnostics(u):
    """(dim Z(A), bip(A), face dimension, number of tight edges outside A), every anti-diagonal covered by A."""
    k = len(u); _, tight = GR.tight_graph(u); A = GR.allowed_edges(k, tight)
    d = GR.dimensions(k, A)
    bip = sum(c["bipartite"] for c in GR.components(k, A))
    assert d["ydim"] == bip and {i + j for i, j in A} == set(range(2 * k - 1))
    return d["zdim"], bip, d["face_dim"], len(tight - A)


def selector_counterexample():
    """The tail (4,5,4,1,5,5,5,5,1,1): K* = {4,5,6}, A* = {6,7,8}, B* empty; block 6 = min A* has no rigid partner,
    while the unique pair {7, 8} at level 15 certifies E-flat."""
    tail = [4, 5, 4, 1, 5, 5, 5, 5, 1, 1]; r = len(tail)
    E, A, pairs, M = attainment(tail)
    hall, K = balanced_by_subsets(A)
    AK = set().union(*(A[q] for q in K)); B = {d for d in range(1, r + 1) if tail[d - 1] == 1 and d - 1 in K}
    rigid = [(s, p[0]) for s, p in enumerate(pairs) if len(p) == 1]
    return (hall and K == {4, 5, 6} and AK == {6, 7, 8} and not B and min(AK) == 6
            and not any(6 in pr for _, pr in rigid) and pairs[15] == [(7, 8)])


def outside_range_example():
    """g = (2, 2,3,2,3,2,3,2,3,2,3, 4): k = 31, r = 11, tail sum 29, superadditive prefix, a tight permutation but
    not total support, and the anti-diagonal 2 E_0 carries no allowed edge."""
    g = [2] + [2, 3] * 5 + [4]
    u = staircase_u(g); k = len(u)
    _, tight = GR.tight_graph(u); A = GR.allowed_edges(k, tight)
    ends, nbs, att, K = tail_structure(g)
    return (k == 31 and len(g) - 1 == 11 and sum(g[1:]) == 29 and A is not None and A != tight
            and superadditive_prefix(g) and bool(K) and not any(i + j == 2 * ends[0] for i, j in A))


def compositions(k):
    for cuts in range(1, (1 << (k - 1)) - 1):
        g, last = [], 0
        for pos in range(1, k):
            if cuts >> (pos - 1) & 1:
                g.append(pos - last); last = pos
        g.append(k - last)
        yield g

