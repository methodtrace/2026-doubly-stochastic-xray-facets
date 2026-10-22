"""Tight and allowed graphs of potentials by forced-edge matching, and exact ranks: a second implementation.

Independent of facet_test.py and facet_criterion.py (which use the switching digraph and modular ranks), for
cross-checks of the newer results. A potential u on {0..k-1} (integers or Fractions):
  tight_graph(u)        (v, E): v(t) = min_{i+j=t} u(i) + u(j), E the ordered tight pairs
  matching(k, E, force) a perfect matching of E (list row -> column) containing the edge force, or None
  allowed_edges(k, E)   the edges lying in some perfect matching (each edge forced in turn), None if there is none
  dimensions(k, G)      zdim = dim Z(G), ydim = dim Y(G) (= number of bipartite components), and the dimension of the
                        face spanned by the X-rays of the perfect matchings of G (when every edge of G is allowed),
                        all by exact rational elimination
  components(k, G)      connected components with a bipartiteness flag
  bipartite_connected   the row/column bipartite graph of E is connected
"""
from fractions import Fraction


def tight_graph(u):
    k = len(u)
    v = [min(u[i] + u[t - i] for i in range(max(0, t - k + 1), min(k - 1, t) + 1)) for t in range(2 * k - 1)]
    return v, {(i, j) for i in range(k) for j in range(k) if u[i] + u[j] == v[i + j]}


def matching(k, edges, force=None):
    adj = [[j for j in range(k) if (i, j) in edges] for i in range(k)]
    xr, xc = force if force is not None else (-1, -1)
    colrow = {}

    def augment(i, seen):
        for j in adj[i]:
            if j == xc or j in seen:
                continue
            seen.add(j)
            if j not in colrow or augment(colrow[j], seen):
                colrow[j] = i
                return True
        return False
    for i in range(k):
        if i != xr and not augment(i, set()):
            return None
    out = [-1] * k
    for j, i in colrow.items():
        out[i] = j
    if force is not None:
        out[xr] = xc
    assert sorted(out) == list(range(k)) and all((i, out[i]) in edges for i in range(k))
    return out


def allowed_edges(k, edges):
    if matching(k, edges) is None:
        return None
    return {e for e in edges if matching(k, edges, e) is not None}


def rank(rows, width):
    piv = {}
    for row in rows:
        cur = [Fraction(x) for x in row]
        assert len(cur) == width
        for c in sorted(piv):
            if cur[c]:
                f = cur[c]
                cur = [a - f * b for a, b in zip(cur, piv[c])]
        p = next((j for j, x in enumerate(cur) if x), None)
        if p is not None:
            d = cur[p]
            piv[p] = [x / d for x in cur]
    return len(piv)


def z_rows(k, edges):
    by = {}
    for i, j in sorted(edges):
        by.setdefault(i + j, []).append((i, j))
    rows = []
    for grp in by.values():
        a, b = grp[0]
        for i, j in grp[1:]:
            row = [0] * k; row[i] += 1; row[j] += 1; row[a] -= 1; row[b] -= 1
            rows.append(row)
    return rows


def zdim(k, edges):
    return k - rank(z_rows(k, edges), k)


def dimensions(k, edges, face=True):
    yrows = []
    for i, j in sorted(edges):
        row = [0] * k; row[i] += 1; row[j] += 1; yrows.append(row)
    out = {"zdim": zdim(k, edges), "ydim": k - rank(yrows, k)}
    if face:
        lst = sorted(edges)
        inc = [[int(i == x) for i, j in lst] for x in range(k)] + [[int(j == x) for i, j in lst] for x in range(k)]
        xr = [[int(i + j == t) for i, j in lst] for t in range(2 * k - 1)]
        out["face_dim"] = rank(inc + xr, len(lst)) - rank(inc, len(lst))
    return out


def components(k, edges):
    adj = [set() for _ in range(k)]
    for i, j in edges:
        adj[i].add(j); adj[j].add(i)
    col, out = {}, []
    for root in range(k):
        if root in col:
            continue
        col[root] = 0; st, verts, bip = [root], [], True
        while st:
            i = st.pop(); verts.append(i)
            for j in adj[i]:
                if j not in col:
                    col[j] = 1 - col[i]; st.append(j)
                elif col[j] == col[i]:
                    bip = False
        out.append({"vertices": sorted(verts), "bipartite": bip})
    return out


def bipartite_connected(k, edges):
    adj = [set() for _ in range(2 * k)]
    for i, j in edges:
        adj[i].add(k + j); adj[k + j].add(i)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x] - seen:
            seen.add(y); st.append(y)
    return len(seen) == 2 * k


def from_increments(d):
    u = [0]
    for x in d:
        u.append(u[-1] + x)
    return u
