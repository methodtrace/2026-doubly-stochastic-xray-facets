"""Grid refinement of facet rays (sec:refinement): spr_m and ref_m, and the lemmas behind the refinement theorem.

For a normalized potential u on {0..K} (k = K + 1) and m >= 2:
    spr_m u (i) = u(ceil(i/m)) on {0..mK},      ref_m u (i) = u(floor(i/m)) on {0..mk - 1},
so every increment of u is followed by m - 1 zeros in spr_m u, and ref_m u = padL^(m-1) spr_m u (the step function of
u read on the grid (1/(mk))Z). With a_i = ceil(i/m) and r_i = m a_i - i:
  Lemma (tight graph of spr_m u): v'(t) = v(ceil(t/m)); (i, j) is tight for spr_m u iff (a_i, a_j) is in E(u) and
      either r_i + r_j <= m - 1, or r_i + r_j >= m and v(a_i + a_j - 1) = v(a_i + a_j).
  Lemma (the graph Theta): for a ray u with connected tight graph and total support, the graph on {1..K} in which
      each set P_s = {p >= 1 : (p, s - p) in E(u)} is a clique is connected.
  Lemma (a flat level): if u is normalized with a tight permutation fixing 0, some (a, b) in E(u) with a, b >= 1
      has v(a + b - 1) = v(a + b).
  Theorem (grid refinement): if u is a normalized ray with connected tight graph and total support, spr_m u is again
      such a ray, gives a facet of hat-P(mK + 1), satisfies P_L, and satisfies P_R if u does; ref_m u gives a facet of
      hat-P(mk).
"""
from .facet_test import tight_edges, v_of


def ceil_div(a, m):
    return -((-a) // m)


def spr(u, m):
    K = len(u) - 1
    return [u[ceil_div(i, m)] for i in range(m * K + 1)]


def ref(u, m):
    return [u[i // m] for i in range(m * len(u))]


def predicted_tight_graph(u, m):
    """The tight graph of spr_m u as predicted by the lemma, from E(u) and v alone."""
    K = len(u) - 1; N = m * K; v = v_of(u); E = tight_edges(u)
    a = [ceil_div(i, m) for i in range(N + 1)]; r = [m * a[i] - i for i in range(N + 1)]
    out = set()
    for i in range(N + 1):
        for j in range(N + 1):
            if (a[i], a[j]) not in E:
                continue
            if r[i] + r[j] <= m - 1 or v[a[i] + a[j] - 1] == v[a[i] + a[j]]:
                out.add((i, j))
    return out


def tight_graph_lemma_holds(u, m):
    K = len(u) - 1; N = m * K; v = v_of(u); up = spr(u, m); vp = v_of(up)
    return all(vp[t] == v[ceil_div(t, m)] for t in range(2 * N + 1)) and tight_edges(up) == predicted_tight_graph(u, m)


def theta_connected(u):
    K = len(u) - 1; E = tight_edges(u); par = list(range(K + 1))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for s in range(2 * K + 1):
        P = sorted({p for p in range(1, K + 1) if 0 <= s - p <= K and (p, s - p) in E})
        for p in P[1:]:
            par[f(p)] = f(P[0])
    return len({f(p) for p in range(1, K + 1)}) == 1


def flat_level(u):
    """Some (a, b) in E(u) with a, b >= 1 and v(a + b - 1) = v(a + b), or None."""
    E = tight_edges(u); v = v_of(u)
    for (a, b) in sorted(E):
        if a >= 1 and b >= 1 and v[a + b - 1] == v[a + b]:
            return (a, b)
    return None


def refined_tight_permutation(sigma, m):
    """The permutation pi of the proof: pi(0) = 0, pi(ma - r) = m sigma(a) - (m - 1 - r)."""
    K = len(sigma) - 1; pi = [0] * (m * K + 1)
    for a in range(1, K + 1):
        for r in range(m):
            pi[m * a - r] = m * sigma[a] - (m - 1 - r)
    return pi
