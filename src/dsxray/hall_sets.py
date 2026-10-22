"""Hall sets of block starts (thm:H).

thm:H: the tight graph of a staircase u(g), 1 <= r <= k-2, fails to have total support iff some nonempty set X
of block starts in {S_1, ..., S_r} has at most |X| tight neighbors.

The condition is decided without enumerating sets: in a bipartite graph, some nonempty row set X has |N(X)| <= |X|
iff there is no matching saturating the rows, or there is none after deleting a single column t.

This module checks the statement on all staircases of given sizes (matching criterion against the certified
total-support test, and against brute force over all subsets for small k), counts minimal Hall sets, and checks
every inequality used in the proof on the actual tight graph.
"""
import itertools
from collections import Counter

from .facet_test import allowed_edges, perfect_matching, tight_edges
from .staircases import ends, staircases, u_of


def saturating(rows, nb, forbid=None):
    """Is there a matching saturating `rows` in the bipartite graph nb (row -> neighbor list), avoiding column forbid?"""
    match = {}

    def aug(s, seen):
        for j in nb[s]:
            if j == forbid or j in seen:
                continue
            seen.add(j)
            if j not in match or aug(match[j], seen):
                match[j] = s
                return True
        return False
    return all(aug(s, set()) for s in rows)


def start_neighbors(g, T=None):
    k = sum(g); r = len(g) - 1; T = T if T is not None else tight_edges(u_of(g))
    E = ends(g); S = [e + 1 for e in E[:-1]]
    return S, {s: sorted(j for j in range(k) if (s, j) in T) for s in S}


def hall_witness_exists(g, T=None):
    """Matching criterion: some nonempty X subset of {S_1..S_r} with |N(X)| <= |X|."""
    k = sum(g); S, nb = start_neighbors(g, T)
    return (not saturating(S, nb)) or any(not saturating(S, nb, t) for t in range(k))


def hall_witness_brute(g, T=None):
    S, nb = start_neighbors(g, T); r = len(S)
    nbs = {s: set(v) for s, v in nb.items()}
    return any(len(set().union(*(nbs[s] for s in X))) <= len(X)
               for n in range(1, r + 1) for X in itertools.combinations(S, n))


def minimal_hall_set_size(g, T=None):
    S, nb = start_neighbors(g, T); nbs = {s: set(v) for s, v in nb.items()}
    for size in range(1, len(S) + 1):
        for X in itertools.combinations(S, size):
            if len(set().union(*(nbs[s] for s in X))) <= size:
                return size
    return None


def check_theorem(g):
    """(total support, matching-criterion witness, neighbors of starts are block ends)."""
    k = sum(g); T = tight_edges(u_of(g)); E = set(ends(g))
    S, nb = start_neighbors(g, T)
    ends_only = all(set(nb[s]) <= E for s in S)
    m = perfect_matching(k, T)
    ts = m is not None and len(allowed_edges(k, T, m)) == len(T)
    return ts, hall_witness_exists(g, T), ends_only


def theorem_counts(k, brute=False):
    """Per-k counts: staircases, lacking total support, exceptions to thm:H, and (optionally) disagreements
    between the matching criterion and brute force over all subsets of starts."""
    C = Counter()
    for g in staircases(k):
        ts, w, eo = check_theorem(g)
        C["staircases"] += 1; C["no_total_support"] += (not ts)
        C["exceptions"] += (ts == w); C["neighbors_not_ends"] += (not eo)
        if brute:
            C["brute_force_disagreements"] += (hall_witness_brute(g) != w)
    return dict(C)


def minimal_size_counts(kmin, kmax):
    C = Counter(); n = 0
    for k in range(kmin, kmax + 1):
        for g in staircases(k):
            u = u_of(g); T = tight_edges(u); m = perfect_matching(k, T)
            if m is not None and len(allowed_edges(k, T, m)) == len(T):
                continue
            n += 1; C[minimal_hall_set_size(g, T)] += 1
    return n, dict(sorted((str(a), b) for a, b in C.items()))


def proof_step_checks(g):
    """Every inequality of the proof of thm:H on the actual tight graph, for every down-closed row set Y given by
    y = (0, y_1..y_r), 0 <= y_p <= g_p:
      (L4) N(S_p) = {E_l : l in A_p}, A_p = {l : m(p, l) = g_p} (lem:starts);
      (L3) neighborhoods of the rows of one block are nested, increasing with the index (fact (F3));
      (L1) the column E_q - (h-1) is a tight neighbor of the row S_p + h - 1 whenever p in P_h, q in A_p, h <= g_q,
           and g_q >= g_p for q >= 1 (Step 4 and fact (F3));
      (L2) |N(Y)| >= sum_h (|A(P_h)| - [h > g_0][0 in A(P_h)]) (the count in Step 4).
    Returns the number of row sets checked (an AssertionError signals a failure)."""
    k = sum(g); r = len(g) - 1; T = tight_edges(u_of(g))
    adj = [set(j for j in range(k) if (i, j) in T) for i in range(k)]
    E = ends(g); S = [0] + [e + 1 for e in E[:-1]]
    M = [max(E[p] + E[s - p] for p in range(max(0, s - r), min(s, r) + 1)) for s in range(2 * r + 1)]
    m = lambda p, q: max(0, E[p] + E[q] - M[p + q - 1])
    A = [{l for l in range(r + 1) if m(p, l) == g[p]} if p >= 1 else None for p in range(r + 1)]
    for p in range(1, r + 1):
        assert adj[S[p]] == {E[l] for l in A[p]}, ("L4", g, p)
        for t in range(g[p] - 1):
            assert adj[S[p] + t] <= adj[S[p] + t + 1], ("L3", g, p)
    n = 0
    for y in itertools.product(*[range(x + 1) for x in g[1:]]):
        y = (0,) + y
        if not any(y):
            continue
        Y = [S[p] + t for p in range(1, r + 1) for t in range(y[p])]
        NY = set().union(*[adj[i] for i in Y]); bound = 0
        for h in range(1, max(y) + 1):
            Ph = [p for p in range(1, r + 1) if y[p] >= h]
            AP = set().union(*[A[p] for p in Ph])
            for p in Ph:
                for q in A[p]:
                    if h <= g[q]:
                        assert E[q] - (h - 1) in adj[S[p] + h - 1], ("L1", g, y)
                    if q >= 1:
                        assert g[q] >= g[p], ("L1", g, y)
            bound += len(AP) - (1 if (h > g[0] and 0 in AP) else 0)
        assert len(NY) >= bound, ("L2", g, y)
        n += 1
    return n
