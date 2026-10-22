"""Exact, output-sensitive enumeration of the rays of the normal fan of hat-P(k): the facet walk.

Mathematics. Q_t = conv{e_i + e_j : i + j = t, i <= j} is a simplex, and the fan of the potential inequalities is the
  (inner) normal fan of the Minkowski sum Q = sum_t Q_t, dim Q = k - 2 (lineality: constants and linear functions).
  Rays of the fan are the facet normals of Q. We walk the facet-ridge graph of Q: from a facet u, the ridges are the
  facets of F(u) = sum_t conv A_t(u), where A_t(u) is the set of pairs attaining v(t); each ridge is crossed by an
  exact line search.

  Facets of a Minkowski sum F of simplices. Let M be the vector matroid of all differences x - x' of vertices lying
  in a common summand; span M is the linear space of F, of dimension D. A facet of F with normal y spans a
  hyperplane H = y^perp of span M, spanned by the differences inside the y-minimal classes; so H is a flat of rank
  D - 1 of M. Conversely, for a rank-(D-1) flat H and each of the two normals y with y^perp = H, y gives a facet iff
  the differences inside the y-minimal classes of the summands span H. The hyperplanes of M are enumerated per
  connected component (the hyperplanes of a direct sum N + N' are H + N' and N + H', H, H' hyperplanes of N, N'), each by a canonical
  depth-first search over flats (a flat is visited only from its lexicographically greedy basis, so each flat once).
  With walk_rank set, a face with at least walk_rank distinct difference directions is handled by walking its own
  facet-ridge graph instead (the recursive variant).

Arithmetic. Ranks are computed modulo p = 2^30 - 35 = 1073741789. Every row is a difference e_i + e_j - e_i' - e_j'
  of two vertices of one simplex Q_t, of squared Euclidean norm at most 6 (the value for i = j) and orthogonal to
  1 and iota, so every rank is at most k - 2, a nonzero minor of that size certifies the rational rank, and by
  Hadamard's inequality every minor of size at most k - 2 is at most 6^((k-2)/2) in absolute value, which is less
  than p for k <= 25; hence the rank mod p equals the rational rank (hadamard_ok asserts 6^(k-2) < p^2). Normals (nullspaces) and line searches use exact Python integers. The facet test in facet_test.py uses a
  different prime, 2^31 - 1, as a one-sided certificate, which needs no bound.

Completeness (prop:complete): the facet graph of a polytope is connected and the ridge
  listing above returns every facet of F(u), so the walk visits every facet of Q. The symmetric-adjacency check
  detects a ridge missed from one side only; it is a consistency check, not a completeness certificate. The
  independent recomputation is the vertex walk (vertex_walk.py, vertex_walk_cdd.py), which shares no code with this
  module.
"""
import sys, random, itertools
from fractions import Fraction
from math import gcd

P = 1073741789          # 2^30 - 35, prime

def vec(x, k):
    r = [0] * k; r[x[0]] += 1; r[x[1]] += 1; return r

def dvec(x, y, k):
    r = [0] * k; r[x[0]] += 1; r[x[1]] += 1; r[y[0]] -= 1; r[y[1]] -= 1; return r

def iprimitive(v):
    """Divide an integer vector by the gcd of its entries."""
    g = 0
    for x in v: g = gcd(g, x)
    return tuple(x // g for x in v) if g > 1 else tuple(v)

def primitive(v):
    den = 1
    for x in v:
        x = Fraction(x); den = den * x.denominator // gcd(den, x.denominator)
    return iprimitive([int(Fraction(x) * den) for x in v])

class Echelon:
    """Incremental row echelon basis mod P (rows normalized to pivot 1)."""
    __slots__ = ("rows",)
    def __init__(self, rows=None): self.rows = list(rows) if rows else []
    def reduce(self, v):
        v = [a % P for a in v]
        for p, r in self.rows:
            c = v[p]
            if c: v = [(a - c * b) % P for a, b in zip(v, r)]
        return v
    def add(self, v):
        v = self.reduce(v)
        for p, a in enumerate(v):
            if a:
                inv = pow(a, P - 2, P); v = [(b * inv) % P for b in v]
                self.rows.append((p, v)); return True
        return False
    def __len__(self): return len(self.rows)

def rank_mod(rows):
    E = Echelon()
    for r in rows: E.add(r)
    return len(E)

class IntEchelon:
    """Fraction-free reduced row echelon form over the integers (each row divided by its content)."""
    __slots__ = ("rows", "k")
    def __init__(self, k, rows=None): self.k = k; self.rows = list(rows) if rows else []   # (pivot, row)
    def copy(self): return IntEchelon(self.k, self.rows)
    def add(self, v):
        v = list(v)
        for p, r in self.rows:
            c = v[p]
            if c:
                a = r[p]; v = [a * x - c * y for x, y in zip(v, r)]
        q = next((i for i, x in enumerate(v) if x), None)
        if q is None: return False
        v = list(iprimitive(v))
        if v[q] < 0: v = [-x for x in v]
        new = []
        for p, r in self.rows:
            c = r[q]
            if c: r = list(iprimitive([v[q] * x - c * y for x, y in zip(r, v)]))
            new.append((p, r))
        new.append((q, v)); self.rows = new
        return True
    def nullspace(self):
        piv = {p: r for p, r in self.rows}; out = []
        for f in range(self.k):
            if f in piv: continue
            L = 1
            for p, r in self.rows: L = L * r[p] // gcd(L, r[p])
            v = [0] * self.k; v[f] = L
            for p, r in self.rows: v[p] = -(L // r[p]) * r[f]
            out.append(iprimitive(v))
        return out

def val(w, x): return w[x[0]] + w[x[1]]

def min_classes(summands, y):
    out = []
    for s in summands:
        vals = [val(y, x) for x in s]; m = min(vals)
        out.append(tuple(x for x, v in zip(s, vals) if v == m))
    return out

def face_rank(summands, k):
    E = Echelon()
    for s in summands:
        for x in s[1:]: E.add(dvec(x, s[0], k))
    return len(E)

def hadamard_ok(k):
    """Minors of size <= k - 2 of rows of squared norm <= 6 are at most 6^((k-2)/2); need that < P (k <= 25)."""
    return 6 ** (k - 2) < P * P

# ---------------- facets of a Minkowski sum of simplices ----------------
class Facets:
    def __init__(self, k, walk_rank=None, seed=4321):
        """walk_rank: a Minkowski sum with at least walk_rank distinct difference directions is handled by walking its
        own facet-ridge graph (recursively, output-sensitive) instead of enumerating the hyperplanes of its matroid
        components (whose number grows like a binomial coefficient in the number of directions); None = never."""
        self.k = k; self.memo = {}; self.walk_rank = walk_rank; self.rng = random.Random(seed)
        self.stats = {"calls": 0, "flats": 0, "hyper": 0, "facets": 0, "walks": 0, "cross": 0}

    def __call__(self, summands):
        S = tuple(sorted(tuple(sorted(s)) for s in summands if len(s) >= 2))
        r = self.memo.get(S)
        if r is not None: return r
        self.stats["calls"] += 1
        res = self._facets(S)
        self.memo[S] = res
        return res

    def _elements(self, S):
        """Distinct difference directions (up to sign/scale)."""
        k = self.k; seen = set(); els = []
        for s in S:
            for x, y in itertools.combinations(s, 2):
                pv = iprimitive(dvec(y, x, k))
                if pv[next(i for i, a in enumerate(pv) if a)] < 0: pv = tuple(-a for a in pv)
                if pv in seen: continue
                seen.add(pv); els.append(list(pv))
        return els

    def _components(self, els):
        """Connected components of the vector matroid (fundamental circuits of a greedy basis), mod P."""
        n = len(els); rows = []; parent = list(range(n)); rk = 0
        def find(a):
            while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for i, e in enumerate(els):
            v = [a % P for a in e]; combo = {}
            for p, r, cb in rows:
                c = v[p]
                if c:
                    v = [(a - c * b) % P for a, b in zip(v, r)]
                    for j, cf in cb.items(): combo[j] = (combo.get(j, 0) - c * cf) % P
            piv = next((p for p, a in enumerate(v) if a), None)
            if piv is None:
                for j, cf in combo.items():
                    if cf: parent[find(j)] = find(i)
            else:
                inv = pow(v[piv], P - 2, P); v = [(a * inv) % P for a in v]
                combo[i] = 1
                combo = {j: (cf * inv) % P for j, cf in combo.items() if cf % P}
                rows.append((piv, v, combo)); rk += 1
        comps = {}
        for i in range(n): comps.setdefault(find(i), []).append(i)
        return [[els[i] for i in c] for c in comps.values()], rk

    def _hyperplanes(self, els, r):
        """All rank-(r-1) flats of the vector matroid els (rank r), each as a list of r-1 independent vectors.
        Canonical DFS: a flat is reached only from its greedy basis; the residuals of the elements outside the
        current flat (mod P) are carried along, so a child costs one elimination step per element."""
        n = len(els)
        if r == 1: return [[]]
        out = []
        root = (tuple(), 0, {f: [a % P for a in els[f]] for f in range(n)})
        stack = [root]; full = (1 << n) - 1
        while stack:
            I, Fmask, res = stack.pop()
            self.stats["flats"] += 1
            last = I[-1] if I else -1
            for e in range(last + 1, n):
                if Fmask >> e & 1: continue
                re = res[e]; q = next(i for i, a in enumerate(re) if a); inv = pow(re[q], P - 2, P)
                re = [(a * inv) % P for a in re]
                F2 = Fmask | (1 << e); nres = {}; ok = True
                for f, rf in res.items():
                    if f == e: continue
                    c = rf[q]
                    if c:
                        rf = [(a - c * b) % P for a, b in zip(rf, re)]
                        if not any(rf):
                            if f < e: ok = False; break
                            F2 |= 1 << f; continue
                    nres[f] = rf
                if not ok: continue
                if len(I) + 1 == r - 1: out.append([els[i] for i in I + (e,)]); self.stats["flats"] += 1
                elif F2 != full: stack.append((I + (e,), F2, nres))
        return out

    def _facets(self, S):
        k = self.k
        els = self._elements(S)
        if not els: return []
        comps, D = self._components(els)
        if self.walk_rank is not None and D >= 2 and len(els) >= self.walk_rank:
            return self._walk(S, D, els)
        res = []
        for ci, comp in enumerate(comps):
            base = IntEchelon(k)
            for cj, c in enumerate(comps):
                if cj != ci:
                    for e in c: base.add(e)
            rc = rank_mod(comp)
            for Hb in self._hyperplanes(comp, rc):
                self.stats["hyper"] += 1
                E = base.copy()
                for h in Hb: E.add(h)
                y = None
                for nv in E.nullspace():
                    if any(sum(a * b for a, b in zip(nv, e)) != 0 for e in comp): y = nv; break
                if y is None: raise RuntimeError("no normal")
                for sg in (1, -1):
                    ys = tuple(sg * a for a in y)
                    if face_rank(min_classes(S, ys), k) == D - 1:
                        res.append(ys); self.stats["facets"] += 1
        return res

    def _walk(self, S, D, els):
        """Facets of the Minkowski sum S (dimension D) by a depth-first walk of its facet-ridge graph."""
        k = self.k; self.stats["walks"] += 1
        E = Echelon(); basis = []
        for e in els:
            if E.add(e): basis.append(e)
            if len(basis) == D: break
        def key(y): return iprimitive([sum(a * b for a, b in zip(y, e)) for e in basis])
        y0 = initial_facet(S, k, self.rng, D)
        seen = {key(y0): y0}; stack = [y0]
        while stack:
            y = stack.pop()
            for z in self(min_classes(S, y)):
                v = cross_ridge(S, y, z); self.stats["cross"] += 1
                kv = key(v)
                if kv in seen: continue
                if face_rank(min_classes(S, v), k) != D - 1: raise RuntimeError("walk produced a non-facet")
                seen[kv] = v; stack.append(v)
        self.stats["facets"] += len(seen)
        return list(seen.values())

# ---------------- the walk on Q ----------------
def Q_summands(k):
    return tuple(tuple((i, t - i) for i in range(k) if i <= t - i < k) for t in range(2 * k - 1))

def normalize(w, k):
    """Integer representative with w(0) = w(k-1) = 0 (constants and linear functions are the lineality)."""
    w0 = w[0]; d = w[-1] - w0
    return iprimitive([(k - 1) * (w[i] - w0) - d * i for i in range(k)])

def first_crossing(S, base, direction, current):
    """Smallest lam = num/den > 0 at which a vertex outside `current` joins the argmin along base + lam*direction."""
    bn = None; bd = 1
    for s, C in zip(S, current):
        c = C[0]; bc = val(base, c); dc = val(direction, c)
        for x in s:
            if x in C: continue
            dx = val(direction, x)
            if dx < dc:
                num = val(base, x) - bc; den = dc - dx
                if num > 0 and (bn is None or num * bd < bn * den): bn, bd = num, den
    return None if bn is None else (bn, bd)

def cross_ridge(Q, u, y):
    B = min_classes(Q, u)
    C = min_classes(B, y)
    lam = first_crossing(Q, u, y, C)
    if lam is not None:
        n, d = lam; return iprimitive([d * a + n * b for a, b in zip(u, y)])
    if min_classes(Q, y) != C: return iprimitive(y)
    mu = first_crossing(Q, y, [-a for a in u], C)
    if mu is None: raise RuntimeError("no crossing found")
    n, d = mu; return iprimitive([d * b - n * a for a, b in zip(u, y)])

def initial_facet(Q, k, rng, D=None):
    if D is None: D = k - 2
    for attempt in range(100):
        w = [rng.randint(-10**6, 10**6) for _ in range(k)]
        while True:
            G = min_classes(Q, w); r = face_rank(G, k)
            if r == D - 1: return iprimitive(w)
            if r >= D: break
            E = IntEchelon(k)
            for s in G:
                for x in s[1:]: E.add(dvec(x, s[0], k))
            N = E.nullspace()
            coef = [rng.randint(-50, 50) for _ in N]
            y = list(iprimitive([sum(c * v[i] for c, v in zip(coef, N)) for i in range(k)]))
            moved = False
            for sg in (1, -1):
                d = [sg * a for a in y]
                lam = first_crossing(Q, w, d, G)
                if lam is not None:
                    n, dd = lam; w = list(iprimitive([dd * a + n * b for a, b in zip(w, d)])); moved = True; break
            if not moved: break
    raise RuntimeError("initial facet not found")

def fan_rays(k, progress=None, seed=12345, use_symmetry=True, walk_rank=None, check_adjacency=True):
    """All rays of the min-plus fan of P(k), normalized u(0) = u(k-1) = 0, primitive integers.

    Completeness: the facet-ridge graph of a polytope is connected, so the walk reaches every facet of Q provided every
    ridge of every visited facet is found (prop:complete). check_adjacency verifies the walk's own
    consistency: every ridge is shared by exactly two facets, so the neighbor relation must be symmetric. This detects a
    ridge missed from one side only; a ridge missed from both sides leaves the relation symmetric."""
    assert hadamard_ok(k), "rank mod P no longer certified exact"
    Q = Q_summands(k); FF = Facets(k, walk_rank); rng = random.Random(seed)
    u0 = normalize(initial_facet(Q, k, rng), k)
    def mirror(u): return tuple(reversed(u)) if use_symmetry else u
    def rep(u): return min(u, mirror(u))
    # walk on mirror orbits: neighbors of mirror(u) are the mirrors of the neighbors of u
    seen = {u0, mirror(u0)}; stack = [rep(u0)]; ncross = 0; nbrs = {}
    while stack:
        u = stack.pop()
        if u in nbrs: continue
        F = min_classes(Q, u); nb = []
        for y in FF(F):
            v = normalize(cross_ridge(Q, u, y), k); ncross += 1; nb.append(v)
            if v in seen: continue
            if face_rank(min_classes(Q, v), k) != k - 3: raise RuntimeError(f"non-facet {v}")
            seen.add(v); seen.add(mirror(v)); stack.append(rep(v))
            if progress and len(seen) // 2 % progress == 0:
                print(f"  rays {len(seen)} stack {len(stack)} memo {len(FF.memo)} {FF.stats}", file=sys.stderr, flush=True)
        if len(set(nb)) != len(nb): raise RuntimeError(f"two ridges of {u} lead to the same facet")
        nbrs[u] = set(nb)
    st = {"cross": ncross, **{("inner_" + a if a == "cross" else a): b for a, b in FF.stats.items()}, "memo": len(FF.memo),
          "orbit_reps": len(nbrs)}
    if check_adjacency:
        bad = 0; ridges = 0
        for u, nb in nbrs.items():
            for v in nb:
                ridges += 1; w = rep(v)
                back = u if w == v else mirror(u)
                if back not in nbrs[w]: bad += 1
        if bad: raise RuntimeError(f"asymmetric facet adjacency: {bad} one-sided ridges")
        st["adjacency_symmetric"] = True; st["ridge_incidences_reps"] = ridges
    return sorted(seen), st

# ---------------- parallel frontier walk (same output as fan_rays) ----------------
_W = {}
def _init_worker(k, walk_rank):
    _W["k"] = k; _W["Q"] = Q_summands(k); _W["FF"] = Facets(k, walk_rank)

def _neighbors(u):
    k = _W["k"]; Q = _W["Q"]; FF = _W["FF"]; nb = []
    for y in FF(min_classes(Q, u)):
        v = normalize(cross_ridge(Q, u, y), k)
        if face_rank(min_classes(Q, v), k) != k - 3: raise RuntimeError(f"non-facet {v}")
        nb.append(v)
    if len(set(nb)) != len(nb): raise RuntimeError(f"two ridges of {u} lead to the same facet")
    return u, nb

def fan_rays_parallel(k, workers, seed=12345, walk_rank=None, log=None, chunksize=4):
    """Breadth-first version of fan_rays over mirror orbits; the frontier is processed by a process pool. Returns the
    same ray set and the same adjacency check."""
    from multiprocessing import Pool
    import time
    assert hadamard_ok(k), "rank mod P no longer certified exact"
    Q = Q_summands(k); rng = random.Random(seed)
    u0 = normalize(initial_facet(Q, k, rng), k)
    def mirror(u): return tuple(reversed(u))
    def rep(u): return min(u, mirror(u))
    seen = {u0, mirror(u0)}; frontier = [rep(u0)]; nbrs = {}; t0 = time.time(); level = 0
    with Pool(workers, initializer=_init_worker, initargs=(k, walk_rank)) as pool:
        while frontier:
            nxt = []
            for u, nb in pool.imap_unordered(_neighbors, frontier, chunksize=chunksize):
                nbrs[u] = set(nb)
                for v in nb:
                    if v not in seen: seen.add(v); seen.add(mirror(v)); nxt.append(rep(v))
            level += 1
            if log: print(f"  level {level}: processed {len(nbrs)} orbit reps, rays {len(seen)}, next frontier {len(nxt)}, "
                          f"{time.time()-t0:.0f}s", file=log, flush=True)
            frontier = sorted(set(nxt))
    bad = 0; ridges = 0
    for u, nb in nbrs.items():
        for v in nb:
            ridges += 1; w = rep(v)
            if (u if w == v else mirror(u)) not in nbrs[w]: bad += 1
    if bad: raise RuntimeError(f"asymmetric facet adjacency: {bad} one-sided ridges")
    return sorted(seen), {"orbit_reps": len(nbrs), "ridge_incidences_reps": ridges, "adjacency_symmetric": True,
                          "bfs_levels": level, "workers": workers}


# ---------------- resumable parallel walk (used for k = 12) ----------------
def fan_rays_resumable(k, workers, state_path, walk_rank=60, seed=12345, log=None):
    """Same walk as fan_rays_parallel, with every processed orbit representative appended to state_path (JSON lines),
    so that an interrupted run resumes where it stopped. Returns (sorted rays, stats); stats["one_sided_ridges"]
    must be 0."""
    import json
    import os
    import time
    from multiprocessing import Pool
    assert hadamard_ok(k), "rank mod P no longer certified exact"
    mirror = lambda u: tuple(reversed(u)); rep = lambda u: min(u, mirror(u))
    nbrs = {}; u0 = None
    if os.path.exists(state_path):
        for line in open(state_path):
            try:
                rec = json.loads(line)
            except ValueError:
                continue                                     # a line cut off by an interruption
            if "u0" in rec:
                u0 = tuple(rec["u0"])
            else:
                nbrs[tuple(rec["u"])] = [tuple(v) for v in rec["nb"]]
    if u0 is None:
        u0 = normalize(initial_facet(Q_summands(k), k, random.Random(seed)), k)
        with open(state_path, "w") as f:
            f.write(json.dumps({"u0": list(u0), "k": k}) + "\n")
    seen = {u0, mirror(u0)}
    for u, nb in nbrs.items():
        seen.add(u); seen.add(mirror(u))
        for v in nb:
            seen.add(v); seen.add(mirror(v))
    todo = sorted({rep(v) for v in seen} - set(nbrs))
    t0 = time.time(); last = t0
    with open(state_path, "a") as f, Pool(workers, initializer=_init_worker, initargs=(k, walk_rank)) as pool:
        while todo:
            nxt = set()
            for u, nb in pool.imap_unordered(_neighbors, todo, chunksize=1):
                nbrs[u] = nb; f.write(json.dumps({"u": list(u), "nb": [list(v) for v in nb]}) + "\n")
                for v in nb:
                    if v not in seen:
                        seen.add(v); seen.add(mirror(v)); nxt.add(rep(v))
                if log and time.time() - last > 60:
                    f.flush(); last = time.time()
                    print(f"  k={k}: {len(nbrs)} orbit representatives processed, {len(seen)} rays seen, "
                          f"{len(nxt)} queued", file=log, flush=True)
            todo = sorted(nxt - set(nbrs))
    bad = sum(1 for u, nb in nbrs.items() for v in nb if (u if rep(v) == v else mirror(u)) not in set(nbrs[rep(v)]))
    return sorted(seen), {"orbit_reps": len(nbrs), "ridge_incidences_reps": sum(len(nb) for nb in nbrs.values()),
                          "one_sided_ridges": bad}


# ---------------- census by smallest step family ----------------
FAMILY_NAMES = {1: "single", 2: "double", 3: "triple", 4: "four-step", 5: "five-step", 6: "six-step"}


def step_family_index(k):
    """Map normalized ray -> smallest r' = min(r, k-1-r) such that the ray is a signed r-step potential."""
    fam = {}
    for r in range(1, k - 1):
        for J in itertools.combinations(range(1, k), r):
            for s in (1, -1):
                key = normalize([s * sum(i >= j for j in J) for i in range(k)], k)
                fam[key] = min(fam.get(key, 99), min(r, k - 1 - r))
    return fam


def census(k, facet_rays, coordinate_facets):
    """Facets of hat-P(k) by the smallest step family that produces them; 'other' = from no step potential."""
    fam = step_family_index(k)
    out = {"coordinate": coordinate_facets}
    for u in facet_rays:
        m = fam.get(tuple(u))
        nm = FAMILY_NAMES.get(m, f"{m}-step") if m else "other"
        out[nm] = out.get(nm, 0) + 1
    return out
