"""The facet criterion for potential inequalities, empty anti-diagonals and the symmetric dual potential.

Notation as in `facet_test`: u a potential on {0..k-1}, v(t) = min_{i+j=t} u(i) + u(j), E(u) the tight graph, A = E*(u)
the allowed graph (tight edges lying in some tight permutation), and for a symmetric set G of pairs

    Z(G) = {z : z(i) + z(j) depends only on i + j along G},      Y(G) = {y : y(i) + y(j) = 0 for (i, j) in G},

and bip(G) the number of bipartite components of the graph on {0..k-1} with an edge {i, j} for every pair of G (a loop
at i when (i, i) is in G; an isolated vertex is a bipartite component). Then dim Y(G) = bip(G), and the facet
criterion reads: if u has a tight permutation, I(u) defines a facet iff

    (i)  every anti-diagonal carries an allowed edge, and      (ii)  dim Z(A) - bip(A) = 3,

and under (i) the face has dimension (2k - 1) - (dim Z(A) - bip(A)).

The symmetric dual potential: if u has a tight permutation sigma and A meets every anti-diagonal, take a function h
constant on the strong components of the switching digraph D_sigma (arc i -> sigma^{-1}(j) for each tight (i, j)) and
strictly decreasing along arcs between components; for small eps > 0, w(i) = u(i) + eps (h(i) - h(sigma^{-1}(i))) / 2
has the same sum as u and w(i) + w(j) >= v(i + j) with equality exactly on A. `symmetric_dual` builds this w in
exact rational arithmetic; `check_symmetric_dual` verifies every conclusion directly.

All ranks here are exact: integer elimination with gcd reduction (`rank_int`). Ranks modulo a prime are used only
as a first pass where a lower bound suffices.
"""
from fractions import Fraction
from math import gcd

from .facet_test import allowed_edges, face_rows, perfect_matching, tight_edges, v_of, z_rows


# ---------------- exact integer rank ----------------
def rank_int(rows):
    """Rank over Q of an integer matrix, by fraction-free elimination with gcd reduction of every new row."""
    rows = [list(map(int, r)) for r in rows if any(r)]
    if not rows:
        return 0
    piv_rows = []          # list of (pivot column, row) with row[pivot] > 0
    for r in rows:
        r = r[:]
        for c, p in piv_rows:
            if r[c]:
                a, b = p[c], r[c]
                g = gcd(a, b); fa, fb = b // g, a // g
                r = [fb * x - fa * y for x, y in zip(r, p)]
        nz = next((c for c, x in enumerate(r) if x), None)
        if nz is None:
            continue
        g = 0
        for x in r:
            if x:
                g = gcd(g, x)
        if r[nz] < 0:
            g = -g
        r = [x // g for x in r]
        piv_rows.append((nz, r))
    return len(piv_rows)


# ---------------- graphs ----------------
def allowed_graph(u):
    """(E, A, match): the tight graph, the allowed graph (None if no tight permutation) and a tight permutation as a
    dict column -> row."""
    k = len(u); E = tight_edges(u); m = perfect_matching(k, E)
    if m is None:
        return E, None, None
    return E, allowed_edges(k, E, m), m


def graph_z_dim(k, G):
    """dim Z(G), exact."""
    rows = z_rows(k, G)
    return k - rank_int(rows)


def y_dim(k, G):
    """dim Y(G) = dim {y : y(i) + y(j) = 0 on G}, exact."""
    rows = []
    for (i, j) in G:
        if i <= j:
            r = [0] * k; r[i] += 1; r[j] += 1; rows.append(r)
    return k - rank_int(rows)


def components(k, G):
    """Connected components of the graph on {0..k-1} with an edge {i, j} per pair of G (loops allowed), each as
    (sorted vertex list, is_bipartite)."""
    adj = [set() for _ in range(k)]
    for (i, j) in G:
        adj[i].add(j); adj[j].add(i)
    col = [-1] * k; out = []
    for s in range(k):
        if col[s] >= 0:
            continue
        col[s] = 0; st = [s]; comp = [s]; bip = True
        while st:
            x = st.pop()
            for y in adj[x]:
                if col[y] < 0:
                    col[y] = 1 - col[x]; st.append(y); comp.append(y)
                elif col[y] == col[x]:
                    bip = False
        out.append((sorted(comp), bip))
    return out


def bip(k, G):
    return sum(1 for _, b in components(k, G) if b)


def empty_antidiagonals(k, G):
    return [t for t in range(2 * k - 1) if not any((i, t - i) in G for i in range(max(0, t - k + 1), min(t, k - 1) + 1))]


def face_dim_exact(k, A):
    """Dimension of the face of hat-P(k) spanned by the X-rays of the perfect matchings of A, assuming every edge of
    A lies in a perfect matching of A: rank [M; X] - rank M (exact)."""
    M, X = face_rows(k, A)
    return rank_int(M + X) - rank_int(M)


def bipartite_connected(k, E):
    """The bipartite graph with rows and columns 0..k-1 and an edge row i - column j per pair of E is connected."""
    adj = {("r", i): set() for i in range(k)}; adj.update({("c", j): set() for j in range(k)})
    for i, j in E:
        adj[("r", i)].add(("c", j)); adj[("c", j)].add(("r", i))
    seen = {("r", 0)}; st = [("r", 0)]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    return len(seen) == 2 * k


# ---------------- the criterion ----------------
def criterion(u, face_dim=False):
    """Quantities of the facet criterion for u. Keys: tight_permutation, total_support, cover (every anti-diagonal
    carries an allowed edge), n_empty (anti-diagonals without allowed edge), zdim (dim Z(u)), zdim_star (dim Z(A)),
    bip (bip(A)), criterion_facet; with face_dim=True also face_dim (exact) and formula_dim."""
    k = len(u); E, A, m = allowed_graph(u)
    out = {"k": k, "tight_permutation": A is not None, "zdim": graph_z_dim(k, E)}
    if A is None:
        out.update(total_support=False, cover=False, criterion_facet=False)
        return out
    empty = empty_antidiagonals(k, A)
    zs = graph_z_dim(k, A); b = bip(k, A)
    out.update(total_support=len(A) == len(E), lost_edges=len(E) - len(A), cover=not empty, n_empty=len(empty),
               zdim_star=zs, bip=b, criterion_facet=(not empty) and zs - b == 3)
    if face_dim:
        out["face_dim"] = face_dim_exact(k, A)
        out["formula_dim"] = (2 * k - 1) - (zs - b) if not empty else None
    return out


# ---------------- the symmetric dual potential ----------------
def _scc(n, succ):
    index = [-1] * n; low = [0] * n; on = [False] * n; st = []; comp = [-1] * n; c = [0, 0]

    def visit(x):
        work = [(x, iter(succ[x]))]; index[x] = low[x] = c[0]; c[0] += 1; st.append(x); on[x] = True
        while work:
            y, it = work[-1]
            for z in it:
                if index[z] < 0:
                    index[z] = low[z] = c[0]; c[0] += 1; st.append(z); on[z] = True
                    work.append((z, iter(succ[z]))); break
                if on[z]:
                    low[y] = min(low[y], index[z])
            else:
                work.pop()
                if work:
                    low[work[-1][0]] = min(low[work[-1][0]], low[y])
                if low[y] == index[y]:
                    while True:
                        z = st.pop(); on[z] = False; comp[z] = c[1]
                        if z == y:
                            break
                    c[1] += 1
    for x in range(n):
        if index[x] < 0:
            visit(x)
    return comp, c[1]


def symmetric_dual(u):
    """The potential w of the symmetric-dual theorem (exact Fractions), or None if u has no tight permutation or
    some anti-diagonal carries no allowed edge."""
    k = len(u); E, A, m = allowed_graph(u)
    if A is None or empty_antidiagonals(k, A):
        return None
    sinv = dict(m)                      # column j -> row sigma^{-1}(j)
    succ = [sorted({sinv[j] for j in range(k) if (i, j) in E}) for i in range(k)]
    comp, nc = _scc(k, succ)
    # h on components: longest path to a sink in the condensation, so every arc between components decreases h
    csucc = [set() for _ in range(nc)]
    for i in range(k):
        for y in succ[i]:
            if comp[y] != comp[i]:
                csucc[comp[i]].add(comp[y])
    hc = [None] * nc

    def H(c):
        if hc[c] is None:
            hc[c] = 1 + max((H(d) for d in csucc[c]), default=-1)
        return hc[c]
    h = [H(comp[i]) for i in range(k)]
    v = v_of(u)
    slack = [u[i] + u[j] - v[i + j] for i in range(k) for j in range(k) if (i, j) not in E]
    delta = min(slack) if slack else 1
    eps = Fraction(delta, 2 * (max(h) + 1))
    sigma_inv = [sinv[j] for j in range(k)]
    return [Fraction(u[i]) + eps * (h[i] - h[sigma_inv[i]]) / 2 for i in range(k)]


def check_symmetric_dual(u, w=None):
    """Verify the conclusions of the symmetric-dual theorem for u (with w = symmetric_dual(u) by default). Returns a
    dict of booleans; every entry is True when the theorem's conclusions hold."""
    k = len(u); E, A, m = allowed_graph(u)
    w = symmetric_dual(u) if w is None else w
    if w is None:
        return {"applicable": False}
    v = v_of(u)
    sl = {(i, j): w[i] + w[j] - v[i + j] for i in range(k) for j in range(k)}
    vw = [min(w[i] + w[t - i] for i in range(max(0, t - k + 1), min(t, k - 1) + 1)) for t in range(2 * k - 1)]
    Ew = {(i, j) for i in range(k) for j in range(k) if w[i] + w[j] == vw[i + j]}
    mw = perfect_matching(k, Ew)
    Aw = allowed_edges(k, Ew, mw) if mw is not None else set()
    d = [u[i] - w[i] for i in range(k)]
    return {"applicable": True,
            "same_sum": sum(w) == sum(u),
            "dominates_v": all(x >= 0 for x in sl.values()),
            "equality_exactly_on_allowed": {e for e, x in sl.items() if x == 0} == A,
            "same_v": vw == v,
            "tight_graph_is_allowed_graph": Ew == A,
            "w_total_support": mw is not None and Aw == Ew,
            "difference_in_Y": all(d[i] + d[j] == 0 for (i, j) in A)}
