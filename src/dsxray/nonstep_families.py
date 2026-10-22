"""Infinite families of facets that come from no step potential (thm:fam) and their checks,
and the family g_b = G_{1,b} of thm:height (height b on 2b + 2 points).

Names in the code and in the paper: family A = u_{a,b} (thm:fam(a)), W = omega (thm:fam(b)), V+ = nu and V- = -nu
(thm:fam(c)).

Family A: pattern 1^a 0^b 2 0 (a >= 2, b >= 1), k = a + b + 3, i.e. u(i) = min(i, a) on [0, n], n = a + b, and
          u(n+1) = u(n+2) = a + 2. One sign.
Family W: u = e_1 - e_{k-2}; facet in both signs (the second sign is the mirror image).
Family V: V+ : u = e_{k-1} - e_0 - e_2 and V- : u = e_0 + e_2 - e_{k-1}; both signs are facets.
Dip:      u = -(e_a + e_{a+2}); facet iff a not in {2, k-5} (checked, not proved; not used in the paper).

For each family the proofs give the tight graph, a universal vertex and a tight permutation explicitly; the
functions below return those predictions so that they can be compared with the computed tight graph.
"""
from .facet_test import classify, is_ray, has_total_support, tight_edges
from .padding import mirror, pad_left, pad_right, pattern


def _sym(E):
    return E | {(j, i) for (i, j) in E}


def family_A(a, b):
    return (1,) * a + (0,) * b + (2, 0)


def family_A_u(a, b):
    n = a + b
    return [min(i, a) for i in range(n + 1)] + [a + 2, a + 2]


def family_A_tight_graph(a, b):
    n = a + b; E = set()
    E |= {(i, j) for i in range(a + 1) for j in range(a + 1) if i + j <= a}
    E |= {(0, t) for t in range(a + 1, n + 1)}
    E |= {(s, n) for s in range(1, a)}
    E |= {(i, j) for i in range(a, n + 1) for j in range(a, n + 1) if i + j >= n + a}
    E |= {(j, n + 2) for j in range(0, a - 1)}
    E |= {(n, n + 1), (n - 1, n + 2), (n, n + 2), (n + 1, n + 2), (n + 2, n + 2)}
    return _sym(E)


def W_u(k):
    u = [0] * k; u[1] = 1; u[k - 2] = -1
    return u


def Vplus_u(k):
    u = [0] * k; u[0] = -1; u[2] = -1; u[k - 1] = 1
    return u


def Vminus_u(k):
    return [-x for x in Vplus_u(k)]


def dip_u(k, a):
    u = [0] * k; u[a] = -1; u[a + 2] = -1
    return u


def W_tight_graph(k):
    Q = range(2, k - 2)
    E = {(0, j) for j in range(k)} | {(1, 0), (1, k - 2)} | {(i, k - 2) for i in range(k)}
    E |= {(i, j) for i in Q for j in Q if i + j <= k - 3 or i + j == k - 1} | {(k - 1, k - 1)}
    return _sym(E)


def Vplus_tight_graph(k):
    Z0 = [1] + list(range(3, k - 1))
    E = {(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (2, 2)}
    E |= {(0, t) for t in range(5, k - 1)} | {(2, j) for j in range(3, k - 1)} | {(2, k - 1)}
    E |= {(i, j) for i in Z0 for j in Z0 if k + 1 <= i + j <= 2 * k - 4} | {(k - 2, k - 1), (k - 1, k - 1)}
    return _sym(E)


def Vminus_tight_graph(k):
    Z0 = [1] + list(range(3, k - 1))
    E = {(0, 0), (0, 1), (1, 1), (0, 3), (1, 2), (0, k - 1), (1, k - 1), (2, k - 1)}
    E |= {(i, j) for i in Z0 for j in Z0 if 4 <= i + j <= k - 1 or i + j == k + 1}
    E |= {(i, k - 1) for i in range(3, k)}
    return _sym(E)


def W_permutation(k):
    return [1, 0] + [k - 1 - i for i in range(2, k - 2)] + [k - 2, k - 1], k - 2


def Vplus_permutation(k):
    return [1, 0, 2] + [k + 1 - i for i in range(3, k - 1)] + [k - 1], 2


def Vminus_permutation(k):
    return [0, 2, 1] + [k + 1 - i for i in range(3, k - 1)] + [k - 1], k - 1


def permutation_with_universal_vertex(E, s, e):
    """s is a tight permutation fixing e, and e is tight with every index."""
    k = len(s)
    return (sorted(s) == list(range(k)) and s[e] == e and all((i, s[i]) in E for i in range(k))
            and all((e, j) in E for j in range(k)))


def family_A_permutation(a, b):
    """Reverse [0, a]; i -> a + n - i on [a+1, n-1]; fix n; swap n+1, n+2."""
    n = a + b
    s = [a - i for i in range(a + 1)] + [a + n - i for i in range(a + 1, n)] + [n, n + 2, n + 1]
    return s, n


def check_member(u):
    """(facet, ray, total support) by the certified test and exact ray test."""
    c = classify(u)
    return c["facet"], is_ray(u), has_total_support(u)


def proved_patterns(k):
    """Patterns of the proved families at k, closed under mirror."""
    out = set(); n = k - 3
    for a in range(2, n):
        out.add(family_A(a, n - a))
    if k >= 9:
        for u in (W_u(k), [-x for x in W_u(k)], Vplus_u(k), Vminus_u(k)):
            out.add(pattern(u))
    return out | {mirror(p) for p in out}


def checked_patterns(k, sparse_templates=None):
    """Exactly checked (not proved) families: Dip, and optional sparse templates {k: [(support, coefficients)]}."""
    out = set()
    for a in range(k - 2):
        if a in (2, k - 5):
            continue
        out.add(pattern(dip_u(k, a)))
    for S, cs in (sparse_templates or {}).get(str(k), []):
        u = [0] * k
        for a, c in zip(S, cs):
            u[a] += c
        out.add(pattern(u))
    return out | {mirror(p) for p in out}


def padding_closure(gen, k, kmin=6):
    """All padL^a padR^b images at k of the patterns gen(j), kmin <= j <= k."""
    mem = set()
    for j in range(kmin, k + 1):
        for p in gen(j):
            for a in range(k - j + 1):
                q = p
                for _ in range(a):
                    q = pad_left(q)
                for _ in range(k - j - a):
                    q = pad_right(q)
                mem.add(q)
    return mem



# ---------------- the family g_b = G_{1,b} of thm:height ----------------
def height_family_pattern(b):
    """Increments of g_b on {0..2b+1}: (b, 0^(b-1), b, 1^b), b >= 2; height b, k = 2b + 2."""
    return (b,) + (0,) * (b - 1) + (b,) + (1,) * b


def height_family_v(b):
    """v for w = g_b - iota as in the proof: 0 at t = 0; b - t on [1, b]; 2b - t on [b+1, 2b]; b - 1 on
    [2b+1, 3b+1]; 2b - 2 on [3b+2, 4b+2]."""
    out = []
    for t in range(4 * b + 3):
        if t == 0:
            out.append(0)
        elif t <= b:
            out.append(b - t)
        elif t <= 2 * b:
            out.append(2 * b - t)
        elif t <= 3 * b + 1:
            out.append(b - 1)
        else:
            out.append(2 * b - 2)
    return out


def height_family_tight_graph(b):
    """Tight pairs of w = g_b - iota as listed in the proof (with transposes): (0, t) for t in [0, b+1] and t = 2b+1;
    pairs in P x P with index sum >= b + 1 (P = [1, b]); (b, j) for j in Q = [b+1, 2b+1]; pairs in Q x Q with index
    sum >= 3b + 2."""
    P = range(1, b + 1); Q = range(b + 1, 2 * b + 2); E = set()
    E |= {(0, t) for t in list(range(0, b + 2)) + [2 * b + 1]}
    E |= {(i, j) for i in P for j in P if i + j >= b + 1}
    E |= {(b, j) for j in Q}
    E |= {(i, j) for i in Q for j in Q if i + j >= 3 * b + 2}
    return _sym(E)


def height_family_permutation(b):
    """The tight permutation of the proof: exchanges 0 and 1, i -> b+1-i on [2, b-1], fixes b, j -> 3b+2-j on Q."""
    k = 2 * b + 2; s = list(range(k))
    s[0], s[1] = 1, 0
    for i in range(2, b):
        s[i] = b + 1 - i
    for j in range(b + 1, 2 * b + 2):
        s[j] = 3 * b + 2 - j
    return s


def quadratic_height_pattern(a, b):
    """Increments (b^a, 0^(b-1), (b-1)a + 1, 1^b), k = a + 2b + 1 (a = 1 gives g_b). An exact check, not a theorem:
    a facet ray iff 1 <= a <= b + 1 (checked for small b); height b^2 at a = b + 1."""
    return (b,) * a + (0,) * (b - 1) + ((b - 1) * a + 1,) + (1,) * b
