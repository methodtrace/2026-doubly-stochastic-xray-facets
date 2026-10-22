"""Statement E-flat and its computer-assisted proof for r <= 6 steps (thm:Eflat6, app:cells).

Attainment graph: for q in {0..r-1}, A(q) = {l in {0..r} : E_q + E_l = M(q+l)}; A(K) is the union over q in K.
Hall's condition: |A(K)| >= |K| for all K. Balanced sets (|A(K)| = |K|) are closed under union under Hall's condition;
K* is the largest balanced set, A* = A(K*), and B* = {d in 1..r : g_d = 1 and d-1 in K*}.

Statement E-flat: if Hall's condition holds and K* is nonempty, some level s is attained by exactly one unordered
pair {c, d} (a rigid level) with c in A* and d not in B*.

E-flat depends only on the tail (g_1, ..., g_r): attainment is decided by the signs of the linear forms
L = E_a + E_b - E_c - E_d (a + b = c + d) in g_1..g_r, and B* by the set I = {d : g_d = 1}. Every integer tail lies in
one real cell {x : x_i = 1 (i in I), x_i >= 2 (i not in I), sign L_j(x) = sigma_j}. `enumerate_cells(r)` lists all
nonempty cells by depth-first branching on the signs, deciding nonemptiness exactly by Fourier-Motzkin elimination
over the rationals with strict and non-strict inequalities, and evaluates E-flat on each cell from the sign data.
"""
import itertools
import sys
from collections import Counter
from fractions import Fraction as F

sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))


# ---------------- E-flat on an integer tail ----------------
def _largest_balanced(G, r):
    """Returns None if Hall's condition fails, else K* (rows that cannot reach a free column by alternating paths)."""
    mc = {}

    def aug(q, seen):
        for c in G[q]:
            if c in seen:
                continue
            seen.add(c)
            if c not in mc or aug(mc[c], seen):
                mc[c] = q
                return True
        return False
    for q in range(r):
        if not aug(q, set()):
            return None
    reach = set(); changed = True
    while changed:
        changed = False
        for q in range(r):
            if q not in reach and any((c not in mc) or (mc[c] in reach) for c in G[q]):
                reach.add(q); changed = True
    return [q for q in range(r) if q not in reach]


def eflat_tail(t):
    """Verdict of E-flat on the tail t = (g_1..g_r): 'nohall', 'ts' (K* empty), 'ok' or 'FAIL'."""
    g = (1,) + tuple(t); r = len(t); E = [0] * (r + 1)
    for l in range(1, r + 1):
        E[l] = E[l - 1] + g[l]
    M = [0] * (2 * r + 1); cnt = [0] * (2 * r + 1)
    for s in range(2 * r + 1):
        b = -1; n = 0
        for p in range(max(0, s - r), s // 2 + 1):
            v = E[p] + E[s - p]
            if v > b:
                b = v; n = 1
            elif v == b:
                n += 1
        M[s] = b; cnt[s] = n
    G = [[l for l in range(r + 1) if E[q] + E[l] == M[q + l]] for q in range(r)]
    Ks = _largest_balanced(G, r)
    if Ks is None:
        return "nohall"
    if not Ks:
        return "ts"
    A = set().union(*[set(G[q]) for q in Ks]); assert len(A) == len(Ks)
    Xe = {d for d in range(1, r + 1) if g[d] == 1 and (d - 1) in Ks}
    for c in A:
        for d in range(r + 1):
            if cnt[c + d] == 1 and E[c] + E[d] == M[c + d] and d not in Xe:
                return "ok"
    return "FAIL"


def compositions_of(n):
    """All compositions of n (tuples)."""
    for mask in range(1 << (n - 1)):
        t = []; c = 1
        for b in range(n - 1):
            if (mask >> b) & 1:
                t.append(c); c = 1
            else:
                c += 1
        t.append(c)
        yield tuple(t)


def eflat_counts(n):
    """E-flat on every tail with g_1 + ... + g_r = n."""
    C = Counter(); bad = []
    for t in compositions_of(n):
        v = eflat_tail(t); C[v] += 1
        if v == "FAIL" and len(bad) < 5:
            bad.append(list(t))
    return {"n": n, "counts": dict(sorted(C.items())), "failures": bad}


# ---------------- exact Fourier-Motzkin feasibility ----------------
def feasible(cons, nvar):
    """Constraints (a, b, kind): a.x + b > 0, >= 0 or == 0 (kind '>', '>=', '='), rational. Exact decision."""
    cons = list(cons)
    while True:
        eq = next((c for c in cons if c[2] == '=' and any(c[0])), None)
        if eq is None:
            break
        a, b, _ = eq; j = next(i for i, x in enumerate(a) if x)
        new = []
        for c in cons:
            if c is eq:
                continue
            ca, cb, kd = c
            if ca[j]:
                f = F(ca[j], 1) / a[j]
                ca = tuple(x - f * y for x, y in zip(ca, a)); cb = cb - f * b
            new.append((ca, cb, kd))
        cons = new
    for c in cons:
        if c[2] == '=' and not any(c[0]) and c[1] != 0:
            return False
    cons = [c for c in cons if c[2] != '=']
    for j in range(nvar):
        pos = []; neg = []; rest = []
        for a, b, kd in cons:
            (pos if a[j] > 0 else neg if a[j] < 0 else rest).append((a, b, kd))
        new = rest[:]
        for (a1, b1, k1) in pos:
            for (a2, b2, k2) in neg:
                f1 = F(-a2[j]); f2 = F(a1[j])
                new.append((tuple(f1 * x + f2 * y for x, y in zip(a1, a2)), f1 * b1 + f2 * b2,
                            '>' if (k1 == '>' or k2 == '>') else '>='))
        dd = {}
        for a, b, kd in new:
            nz = [abs(x) for x in a if x]
            sc = max(nz) if nz else F(1)
            key = (tuple(F(x) / sc for x in a), F(b) / sc)
            if key in dd:
                if kd == '>':
                    dd[key] = '>'
            else:
                dd[key] = kd
        cons = [(a, b, kd) for (a, b), kd in dd.items()]
        for a, b, kd in cons:
            if not any(a) and ((kd == '>' and not b > 0) or (kd == '>=' and not b >= 0)):
                return False
    return all(not ((kd == '>' and not b > 0) or (kd == '>=' and not b >= 0)) for a, b, kd in cons)


# ---------------- cell enumeration ----------------
class CellEnumeration:
    """All nonempty cells for r steps and the E-flat verdict on each."""

    def __init__(self, r):
        self.r = r
        Evec = lambda a: tuple(1 if 1 <= i <= a else 0 for i in range(1, r + 1))
        self.levelpairs = [[(p, s - p) for p in range(max(0, s - r), s // 2 + 1)] for s in range(2 * r + 1)]
        forms = {}; self.flist = []; self.pairform = {}
        for s in range(2 * r + 1):
            for (a, b), (c, d) in itertools.combinations(self.levelpairs[s], 2):
                v = tuple(x + y - z - w for x, y, z, w in zip(Evec(a), Evec(b), Evec(c), Evec(d)))
                if not any(v):
                    raise RuntimeError("a comparison form vanishes identically")
                sg = next(1 if x > 0 else -1 for x in v if x)
                cv = tuple(sg * x for x in v)
                if cv not in forms:
                    forms[cv] = len(self.flist); self.flist.append(cv)
                self.pairform[((a, b), (c, d))] = (forms[cv], sg)
        self.cells = {}

    def verdict(self, sign, I):
        r = self.r; pairform = self.pairform; levelpairs = self.levelpairs

        def cmp(P, Q):
            if P == Q:
                return 0
            if (P, Q) in pairform:
                i, sg = pairform[(P, Q)]; return sg * sign[i]
            i, sg = pairform[(Q, P)]; return -sg * sign[i]
        att = []; uniq = []
        for s in range(2 * r + 1):
            Ps = levelpairs[s]
            A = [P for P in Ps if all(cmp(P, Q) >= 0 for Q in Ps)]
            att.append(set(A)); uniq.append(len(A) == 1)
        attains = lambda q, l: tuple(sorted((q, l))) in att[q + l]
        G = [[l for l in range(r + 1) if attains(q, l)] for q in range(r)]
        Ks = _largest_balanced(G, r)
        if Ks is None:
            return "nohall"
        if not Ks:
            return "ts"
        A = set().union(*[set(G[q]) for q in Ks]); assert len(A) == len(Ks)
        Xe = {d for d in range(1, r + 1) if d in I and (d - 1) in Ks}
        for c in A:
            for d in range(r + 1):
                if uniq[c + d] and tuple(sorted((c, d))) in att[c + d] and d not in Xe:
                    return "ok"
        return "FAIL"

    def run(self, progress=None, feasibility=None):
        """Enumerate the nonempty cells. feasibility(cons, nvar) defaults to the exact routine `feasible`; a wrapper
        can be passed to cross-check the discarded branches."""
        r = self.r; nf = len(self.flist); feas = feasibility or feasible
        for I in itertools.product((0, 1), repeat=r):
            Iset = frozenset(i + 1 for i in range(r) if I[i])
            free = [i for i in range(r) if not I[i]]; nv = len(free)

            def tocons(form, kind):
                return (tuple(F(form[i]) for i in free), F(sum(form[i] for i in range(r) if I[i])), kind)
            base = [(tuple(F(1) if j == t else F(0) for j in range(nv)), F(-2), '>=') for t in range(nv)]
            sign = [0] * nf

            def dfs(i, cons):
                if i == nf:
                    self.cells[(Iset, tuple(sign))] = self.verdict(sign, Iset)
                    return
                f = self.flist[i]
                for sg, kind, ff in ((1, '>', f), (-1, '>', tuple(-x for x in f)), (0, '=', f)):
                    c = tocons(ff, kind)
                    if not any(c[0]):
                        if not ((kind == '>' and c[1] > 0) or (kind == '=' and c[1] == 0)):
                            continue
                        sign[i] = sg; dfs(i + 1, cons)
                    else:
                        nc = cons + [c]
                        if feas(nc, nv):
                            sign[i] = sg; dfs(i + 1, nc)
            dfs(0, base)
            if progress:
                progress(sorted(Iset), len(self.cells))
        return self.cells

    def summary(self):
        C = Counter(self.cells.values())
        return {"r": self.r, "forms": len(self.flist), "cells": len(self.cells), "hall_fails": C["nohall"],
                "kstar_empty": C["ts"], "eflat_holds": C["ok"], "eflat_fails": C["FAIL"]}

    def cell_of(self, t):
        sg = []
        for f in self.flist:
            v = sum(a * b for a, b in zip(f, t)); sg.append((v > 0) - (v < 0))
        return (frozenset(i + 1 for i in range(self.r) if t[i] == 1), tuple(sg))

    def safeguard(self, nmax):
        """Every integer tail with r parts and sum <= nmax lies in an enumerated cell whose verdict equals the
        verdict computed directly from the integers. Returns (tails, distinct cells realized, missing)."""
        def parts(n, r):
            if r == 1:
                yield (n,); return
            for a in range(1, n - r + 2):
                for rest in parts(n - a, r - 1):
                    yield (a,) + rest
        seen = set(); miss = 0; n_t = 0
        for n in range(self.r, nmax + 1):
            for t in parts(n, self.r):
                n_t += 1; key = self.cell_of(t); seen.add(key)
                if key not in self.cells:
                    miss += 1
                else:
                    assert eflat_tail(t) == self.cells[key], (t, eflat_tail(t), self.cells[key])
        return n_t, len(seen), miss

    def export(self):
        """Cells as JSON-ready records: unit steps I, sign vector, verdict; plus the list of forms."""
        rows = [{"unit_steps": sorted(I), "signs": list(s), "verdict": v} for (I, s), v in sorted(
            self.cells.items(), key=lambda kv: (sorted(kv[0][0]), kv[0][1]))]
        return {"r": self.r, "forms": [list(f) for f in self.flist], "cells": rows}


# ---------------- consistency with the tight graph ----------------
def largest_closed_start_set_matches(g):
    """For a staircase with a tight permutation: K* from the attainment graph equals the largest set of starts closed
    under the switching digraph on the actual tight graph; and E-flat's verdict agrees with total support.
    Returns None (no tight permutation), 'ts' or 'nts'; raises AssertionError on disagreement."""
    from .facet_test import perfect_matching, tight_edges
    from .staircases import ends, u_of
    k = sum(g); r = len(g) - 1; T = tight_edges(u_of(g))
    adj = [[j for j in range(k) if (i, j) in T] for i in range(k)]
    m = perfect_matching(k, T)
    if m is None:
        return None
    rowof = m                                      # column -> row
    out = [[rowof[j] for j in adj[i]] for i in range(k)]
    comp = _scc(k, out)
    ts = all(comp[i] == comp[rowof[j]] for i in range(k) for j in adj[i])
    E = ends(g); starts_ = {E[l] + 1 for l in range(r)}
    nc = max(comp) + 1; cout = [set() for _ in range(nc)]
    for i in range(k):
        for j in adj[i]:
            if comp[i] != comp[rowof[j]]:
                cout[comp[i]].add(comp[rowof[j]])
    mem = [[i for i in range(k) if comp[i] == c] for c in range(nc)]

    def clos(c0):
        S = {c0}; st = [c0]
        while st:
            v = st.pop()
            for w in cout[v]:
                if w not in S:
                    S.add(w); st.append(w)
        return S
    Xg = {x for c in range(nc) if all(set(mem[d]) <= starts_ for d in clos(c)) for x in mem[c]}
    M = [max(E[p] + E[s - p] for p in range(max(0, s - r), min(s, r) + 1)) for s in range(2 * r + 1)]
    G = [[l for l in range(r + 1) if E[q] + E[l] == M[q + l]] for q in range(r)]
    Ks = _largest_balanced(G, r)
    assert Ks is not None, g
    assert Xg == {E[q] + 1 for q in Ks}, g
    v = eflat_tail(g[1:])
    assert (v == "ts") == ts and (ts or v == "ok"), (g, v, ts)
    return "ts" if ts else "nts"


def _scc(n, out):
    index = [-1] * n; low = [0] * n; on = [False] * n; st = []; comp = [-1] * n; counter = [0, 0]
    for s0 in range(n):
        if index[s0] >= 0:
            continue
        stack = [(s0, 0)]; index[s0] = low[s0] = counter[0]; counter[0] += 1; st.append(s0); on[s0] = True
        while stack:
            v, pi = stack[-1]
            if pi < len(out[v]):
                stack[-1] = (v, pi + 1); w = out[v][pi]
                if index[w] < 0:
                    index[w] = low[w] = counter[0]; counter[0] += 1; st.append(w); on[w] = True; stack.append((w, 0))
                elif on[w]:
                    low[v] = min(low[v], index[w])
            else:
                stack.pop()
                if stack:
                    low[stack[-1][0]] = min(low[stack[-1][0]], low[v])
                if low[v] == index[v]:
                    while True:
                        w = st.pop(); on[w] = False; comp[w] = counter[1]
                        if w == v:
                            break
                    counter[1] += 1
    return comp


def lp_feasible_float(cons, nvar):
    """FLOAT cross-check of `feasible` (scipy HiGHS): maximize a slack s <= 1 with a.x + b >= s for strict, a.x + b >= 0
    for non-strict constraints and a.x + b = 0 for equations; feasible iff the optimum is positive (or unbounded)."""
    from scipy.optimize import linprog
    A_ub = []; b_ub = []; A_eq = []; b_eq = []
    for a, b, k in cons:
        a = [float(x) for x in a]; b = float(b)
        if k == '=':
            A_eq.append(a + [0.0]); b_eq.append(-b)
        elif k == '>=':
            A_ub.append([-x for x in a] + [0.0]); b_ub.append(b)
        else:
            A_ub.append([-x for x in a] + [1.0]); b_ub.append(b)
    res = linprog(c=[0.0] * nvar + [-1.0], A_ub=A_ub or None, b_ub=b_ub or None, A_eq=A_eq or None, b_eq=b_eq or None,
                  bounds=[(None, None)] * nvar + [(None, 1.0)], method="highs")
    if res.status == 3:
        return True
    return res.status == 0 and -res.fun > 1e-9


def random_systems_crosscheck(n, seed=0):
    """`feasible` against the floating-point LP on n random small systems (1 to 4 variables, 1 to 6 constraints with
    coefficients in -2..2). Returns (agreements, disagreements)."""
    import random
    rnd = random.Random(seed); agree = dis = 0
    for _ in range(n):
        nv = rnd.randint(1, 4); m = rnd.randint(1, 6); cons = []
        for _ in range(m):
            a = tuple(F(rnd.randint(-2, 2)) for _ in range(nv)); b = F(rnd.randint(-3, 3))
            cons.append((a, b, rnd.choice(['>', '>=', '='])))
        if feasible(cons, nv) == lp_feasible_float(cons, nv):
            agree += 1
        else:
            dis += 1
    return agree, dis
