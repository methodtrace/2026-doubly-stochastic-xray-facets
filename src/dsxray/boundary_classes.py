"""Classes of boundary steps, statement (S1) and lem:R (a solitary index excludes total support).

Notation: E_0 < ... < E_r the block ends of u(g), M(s) = max_{p+q=s} E_p + E_q, P(q) = {l : E_q + E_l = M(q+l)} the
attaining partners of q, and the classes of {0..r-1} (two indices are related when they occur in attaining pairs of
a common level). An index q < r is solitary if |P(q)| <= 1.

lem:R: if some q in {0..r-1} is solitary, the tight graph of u(g) lacks total support.
Statement (S1): if the boundary steps form at least two classes, some q in {0..r-1} is solitary.
(S1) and lem:R give: total support implies one class for g and for g*, hence the ray property (thm:S1).

All functions use integer arithmetic. The C program tools/c/boundary_classes.c implements `s1_holds` for the
exhaustive ranges of the paper; `s1_counts` is the Python fallback for small k.
"""
import itertools
from fractions import Fraction

from .level_graph import attain_matrix, boundary_classes, hplus_verdict
from .staircases import ends, staircases, u_of


def solitary_indices(E, att=None):
    r = len(E) - 1; att = att or attain_matrix(E)
    return [q for q in range(r) if sum(att[q]) <= 1]


def s1_holds(E):
    """Returns (several classes, (S1) holds, a solitary index lies outside the class of 0)."""
    att = attain_matrix(E); cl = boundary_classes(E, att)
    if len(cl) == 1:
        return False, True, False
    sol = solitary_indices(E, att)
    c0 = next(c for c in cl if 0 in c)
    return True, bool(sol), any(q not in c0 for q in sol)


def s1_counts(k):
    n = multi = fail = in_u = 0
    for g in staircases(k):
        n += 1
        several, ok, outside = s1_holds(ends(g))
        if several:
            multi += 1; fail += (not ok); in_u += outside
    return {"k": k, "compositions": n, "several_classes": multi, "solitary_outside_class_of_0": in_u, "s1_failures": fail}


def s1_pattern_lifts(V, r):
    """(S1) on all sequences E_p = P_p + (V+1) p with P in {0..V}^{r+1}, min P = 0 (by affine invariance these are
    end sequences of compositions)."""
    multi = fail = 0
    for P in itertools.product(range(V + 1), repeat=r + 1):
        if min(P) != 0:
            continue
        E = [P[p] + (V + 1) * p for p in range(r + 1)]
        several, ok, _ = s1_holds(E)
        if several:
            multi += 1; fail += (not ok)
    return {"V": V, "r": r, "several_classes": multi, "s1_failures": fail}


# ---------------- direct tight graph (from v = u box u, not from the corner lemma) ----------------
def tight_direct(u):
    k = len(u); v = [min(u[i] + u[t - i] for i in range(max(0, t - k + 1), min(t, k - 1) + 1)) for t in range(2 * k - 1)]
    return [[j for j in range(k) if u[i] + u[j] == v[i + j]] for i in range(k)]


def total_support_direct(adj):
    """(has total support, number of tight edges in no tight permutation) by matching and reachability."""
    import sys
    sys.setrecursionlimit(max(100000, sys.getrecursionlimit()))
    k = len(adj); match = {}

    def aug(i, seen):
        for j in adj[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in match or aug(match[j], seen):
                match[j] = i
                return True
        return False
    for i in range(k):
        if not aug(i, set()):
            return False, None
    succ = [{match[j] for j in adj[i]} for i in range(k)]
    reach = []
    for s in range(k):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in succ[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        reach.append(seen)
    bad = sum(1 for i in range(k) for j in adj[i] if not (match[j] == i or (match[j] in reach[i] and i in reach[match[j]])))
    return bad == 0, bad


def dim_z_direct(adj):
    k = len(adj); byt = {}
    for i in range(k):
        for j in adj[i]:
            if i <= j:
                byt.setdefault(i + j, []).append((i, j))
    rows = []
    for L in byt.values():
        i0, j0 = L[0]
        for (i, j) in L[1:]:
            row = [Fraction(0)] * k; row[i] += 1; row[j] += 1; row[i0] -= 1; row[j0] -= 1; rows.append(row)
    rank = 0
    for c in range(k):
        piv = next((rr for rr in range(rank, len(rows)) if rows[rr][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for rr in range(len(rows)):
            if rr != rank and rows[rr][c] != 0:
                f = rows[rr][c] / rows[rank][c]; rows[rr] = [x - f * y for x, y in zip(rows[rr], rows[rank])]
        rank += 1
    return k - rank


def describe(g):
    """Exact description of one composition: classes, (H+) verdict (not used in the paper), (S1), total support and dim Z from the tight
    graph computed directly from u."""
    E = ends(g); att = attain_matrix(E)
    adj = tight_direct(u_of(g))
    ts, unused = total_support_direct(adj)
    several, s1, _ = s1_holds(E)
    return {"g": list(g), "k": sum(g), "r": len(g) - 1, "ends": E, "classes": boundary_classes(E, att),
            "hplus": hplus_verdict(list(g))[0], "s1_holds": s1, "solitary": solitary_indices(E, att),
            "total_support": ts, "unused_tight_edges": unused,
            "tight_edges": sum(len(a) for a in adj), "dim_z": dim_z_direct(adj)}


def lemma_r_counts(kmin, kmax):
    """lem:R on all compositions: those with a solitary index, and how many of them have total support (must be 0)."""
    n = hit = bad = 0
    for k in range(kmin, kmax + 1):
        for g in staircases(k):
            n += 1; E = ends(g)
            if solitary_indices(E):
                hit += 1
                if total_support_direct(tight_direct(u_of(g)))[0]:
                    bad += 1
    return {"compositions": n, "with_solitary_index": hit, "of_which_total_support": bad}
