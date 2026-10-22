#!/usr/bin/env python3
"""Worked examples and small statements of the paper, each checked exactly.

Preliminaries
  - rem:multiplicity (the tail-sum reading): for random step potentials with repeated steps and random lattice laws,
    I(u)/k equals the tail-sum bound with thresholds j/k of multiplicity d(j) (2,000 random cases).
  - After lem:T (sec:SM-prelim-ex): for k = 7, J = {1, 2, 5}, sign +, the tight graph would give face dimension
    10 = 2k - 4, the allowed graph gives 1, and the inequality is not a facet.
  - lem:coord and sec:SM-prelim-ex: for k = 3 the inequality x(2) >= 0 is not a facet; for 4 <= k <= 12 all 2k - 1
    coordinate inequalities are facets.
Padding and families
  - Before cor:cores (sec:padding): the facet ray with increments (0, 2, 0, 2, 1, 1, 1, 0, 1, 1, 1) (k = 12) is a
    core; removing its first increment gives a potential that is not a ray.
  - prop:P on the complete facet lists: for k = md <= 12 (m, d >= 2) no facet normal has the form
    A + Bt + psi(t mod d); at k = 6, exactly two of the four facets from no step potential have normals
    A + Bt + psi(t mod 5).
Step potentials
  - Symmetries (sec:steps, after lem:starts): at k = 7, (2, 1, 1, 3)* = (1, 1, 4, 1), both give facets, and the facets
    are mirror images of each other but not equal.
  - After prop:Sprime: single steps give 2(k-1) distinct facets (4 <= k <= 12), double steps 2 binom(k-1, 2)
    distinct facets (6 <= k <= 12), and at k = 5 the twelve signed double steps give six distinct facets.
  - prop:hinges: v(t) = (t - 2j)^+ for h_j, the folded formula for -h_j, the second differences delta_{2j} and
    -delta_j + delta_{k-1} - delta_{k-1+j}, and the 2(k-2) hinge facets are pairwise distinct (3 <= k <= 16).
  - After cor:smallatoms and rem:LLYtest: at k = 12, 7,008 of the 7,018 ray facets are not convex, and 6,998 are
    neither convex nor concave (data/facets).
  - sec:SM-square (lem:square(ii)): g = (1, 1, 5, 3, 2, 2, 3, 5) (k = 22) has E = (0, 1, 6, 9, 11, 13, 16, 21), block
    4 attains only with itself and with block 7, the pairs (1, 7), (2, 6), (3, 5), (4, 4) attain M(8) = 22; it has
    total support and gives a facet.
Total support of staircases
  - sec:s_census-dict (minimal balanced sets): g = (1, 2, 2, 1, 2, 1, 2, 3, 1, 1, 3, 2, 1, 2, 1, 2, 2, 4, 2, 2)
    (k = 37): the start S_7 = 11 has the single tight neighbor E_10 = 18 (from the tight graph of u(g)), {S_7} is a
    minimal balanced set, no rigid level is attained by a pair containing 10, the largest balanced set is
    {S_6, S_7, S_10}, and exactly nine rigid levels certify it (attaining pair {c, d} with c in A*, d not in B*).
  - sec:s_census-proofs (after the proof of prop:stairbip): g = (1, 1, 4, 4, 7, 4, 4, 1, 6, 5) (k = 37): the classes
    of boundary steps are {E_0, ..., E_7} and {E_8}, the only solitary index is 7, A(7) = {0}; the tight graph has a
    tight permutation but no total support.
Default and --full are the same (about 20 s).
"""
import argparse
import itertools
from functools import reduce
from math import comb, gcd

from common import Timer, finish, load_facet_rows, load_rays, out_dir, report

from dsxray import staircases as ST
from dsxray.cell_enumeration import _largest_balanced
from dsxray.facet_test import (classify, classify_coordinate, face_rows, has_total_support, matchable_edges,
                               perfect_matching, rank_exact, tight_edges, v_of)
from dsxray.level_graph import Mfun, attain_matrix, attainers, boundary_classes
from dsxray.boundary_classes import solitary_indices
from dsxray.padding import from_pattern, pattern
from dsxray.small_k_fan import line_enumeration


def u_J(J, k):
    return [sum(i >= j for j in J) for i in range(k)]


def d2(v):
    return [v[t + 1] - 2 * v[t] + v[t - 1] for t in range(1, len(v) - 1)]


def normal_key(u):
    d = d2(v_of(list(u))); g = reduce(gcd, [abs(x) for x in d]) or 1
    return tuple(x // g for x in d)


def periodic_affine(v, d):
    """v(t) = A + B t + psi(t mod d) for some A, B, psi  <=>  v(t + d) - v(t) is constant."""
    diffs = {v[t + d] - v[t] for t in range(len(v) - d)}
    return len(diffs) <= 1


def tail_sum_xi(eta):
    """Thresholds xi_1 <= ... <= xi_2r of the tail-sum bound of the staircase with thresholds eta (closed formula)."""
    from fractions import Fraction as Fr
    r = len(eta); E = list(eta) + [Fr(1)]
    return [max(E[i] + E[l - 1 - i] for i in range(0, l) if i <= r and l - 1 - i <= r) for l in range(1, 2 * r + 1)]


def tail_sum_violation(eta, law):
    """sum_l P(S >= xi_l) - 2 sum_q (1 - eta_q) for a finite law [(atom, mass), ...]."""
    X = tail_sum_xi(eta)
    return sum(m * sum(1 for x in X if x <= a) for a, m in law) - 2 * sum(1 - e for e in eta)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/examples")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}

    # ---- tail-sum reading of a step potential with repeated thresholds (rem:multiplicity)
    from fractions import Fraction as Fr
    import random
    rng = random.Random(7); bad = n = 0
    for _ in range(2000):
        k = rng.randint(3, 9); r = rng.randint(1, 6)
        J = sorted(rng.randint(1, k - 1) for _ in range(r))
        u = [sum(1 for j in J if j <= i) for i in range(k)]
        law = [(Fr(tt + 1, k), Fr(rng.randint(0, 4))) for tt in range(2 * k - 1)]
        tot = sum(m for _, m in law)
        if tot == 0:
            continue
        law = [(x, m / tot) for x, m in law]
        v = v_of(u)
        lhs = sum(m * v[tt] for tt, (_, m) in enumerate(law)) - Fr(2, k) * sum(u)
        n += 1; bad += lhs != tail_sum_violation([Fr(j, k) for j in J], law)
    checks["tail_sum_multiplicities"] = report("rem:multiplicity: the tail-sum bound with repeated thresholds equals "
                                               "I(u)/k on the grid", bad == 0, f"{n} random cases")

    # ---- preliminaries, padding and families
    k = 7; u = u_J((1, 2, 5), k); E = sorted(tight_edges(u)); Ep = sorted(matchable_edges(k, set(E)))

    def dim(ed):
        M, A = face_rows(k, ed)
        return rank_exact(M + A) - rank_exact(M)
    res["lemma_5_2_example"] = {"tight_graph_dim": dim(E), "allowed_graph_dim": dim(Ep)}
    checks["lemma_5_2_example"] = report("k = 7, J = {1, 2, 5}: tight graph 10 = 2k - 4, allowed graph 1, not a facet",
                                         dim(E) == 10 and dim(Ep) == 1 and not classify(u)["facet"])
    ok = not classify_coordinate(3, 2)["facet"]
    ok &= all(classify_coordinate(kk, t)["facet"] for kk in range(4, 13) for t in range(2 * kk - 1))
    checks["coordinate_facets"] = report("x(2) >= 0 is not a facet for k = 3; all coordinate inequalities are facets "
                                         "for 4 <= k <= 12", ok)
    p = (0, 2, 0, 2, 1, 1, 1, 0, 1, 1, 1); up = from_pattern(p); c = classify(up)
    w = classify(from_pattern(p[1:]))
    res["core_example"] = {"facet": c["facet"], "normalized": pattern(up) == p, "first_increment_removed_is_ray": w["ray"],
                           "last_increment_is_max": p[-1] == max(p)}
    checks["core_example"] = report("k = 12: increments (0,2,0,2,1,1,1,0,1,1,1) give a facet ray that is a core; "
                                    "removing the first increment gives a non-ray",
                                    c["facet"] and c["ray"] and pattern(up) == p and not w["ray"] and p[-1] != max(p))

    t = Timer(); ok = True; per = {}
    for kk in range(4, 13):
        if kk >= 8:
            rows = [r["u"] for r in load_facet_rows(kk) if r["facet"]]
        else:
            rays = load_rays(kk) if kk >= 5 else line_enumeration(kk)[0]
            rows = [list(x) for x in rays if classify(list(x))["facet"]]
        normals = [v_of(list(x)) for x in rows]
        coords = [[-int(s == t0) for s in range(2 * kk - 1)] for t0 in range(2 * kk - 1)]
        for d in range(2, kk):
            if kk % d == 0:
                hits = sum(periodic_affine(v, d) for v in normals + coords)
                per[f"k={kk}, d={d}"] = hits; ok &= hits == 0
    F6 = [list(x) for x in load_rays(6) if classify(list(x))["facet"]]
    nonstep6 = [x for x in F6 if max(pattern(x)) > 1]
    period5 = sum(periodic_affine(v_of(x), 5) for x in nonstep6)
    res["proposition_5_1"] = {"divisor_periodic_facets": per, "k6_nonstep_facets": len(nonstep6), "k6_period5": period5}
    checks["proposition_5_1"] = report("no divisor-periodic facet normal for k = md <= 12; at k = 6, 2 of the 4 non-step "
                                       "facets have period 5", ok and len(nonstep6) == 4 and period5 == 2, f"({t()} s)")

    # ---- step potentials
    g, gs = (2, 1, 1, 3), (1, 1, 4, 1)
    na, nb = normal_key(ST.u_of(g)), normal_key(ST.u_of(gs))
    ok = ST.star(g) == gs and classify(ST.u_of(g))["facet"] and classify(ST.u_of(gs))["facet"]
    ok &= nb == tuple(reversed(na)) and na != nb
    checks["mirror_example"] = report("k = 7: (2,1,1,3)* = (1,1,4,1), mirror-image facets, not equal", ok)

    ok = True; counts = {}
    for kk in range(4, 13):
        singles = {normal_key([s * x for x in u_J((j,), kk)]) for j in range(1, kk) for s in (1, -1)}
        counts[f"single k={kk}"] = len(singles); ok &= len(singles) == 2 * (kk - 1)
    for kk in range(5, 13):
        dbl = [(s, J) for J in itertools.combinations(range(1, kk), 2) for s in (1, -1)]
        keys = {normal_key([s * x for x in u_J(J, kk)]) for s, J in dbl}
        counts[f"double k={kk}"] = len(keys)
        ok &= len(keys) == (6 if kk == 5 else 2 * comb(kk - 1, 2))
        ok &= kk != 5 or len(dbl) == 12
    res["step_distinctness"] = counts
    checks["step_distinctness"] = report("single steps 2(k-1) distinct facets (k >= 4), double steps 2 binom(k-1,2) "
                                         "(k >= 6), six at k = 5", ok)

    ok = True
    for kk in range(3, 17):
        keys = set()
        for j in range(1, kk - 1):
            h = [max(0, i - j) for i in range(kk)]
            ok &= v_of(h) == [max(0, t - 2 * j) for t in range(2 * kk - 1)]
            ok &= d2(v_of(h)) == [int(t == 2 * j) for t in range(1, 2 * kk - 2)]
            fold = [-max(t - j, 0) if t <= kk - 1 else -(kk - 1 - j) - max(t - kk + 1 - j, 0) for t in range(2 * kk - 1)]
            ok &= v_of([-x for x in h]) == fold
            ok &= d2(fold) == [-int(t == j) + int(t == kk - 1) - int(t == kk - 1 + j) for t in range(1, 2 * kk - 2)]
            keys |= {normal_key(h), normal_key([-x for x in h])}
            if kk <= 12:
                ok &= classify(h)["facet"] and classify([-x for x in h])["facet"]
        ok &= len(keys) == 2 * (kk - 2)
    checks["hinges"] = report("hinges: v formulas, second differences, 2(k-2) distinct facets (3 <= k <= 16)", ok)

    t = Timer(); F12 = [r["u"] for r in load_facet_rows(12) if r["facet"]]
    inc = lambda x: [x[i + 1] - x[i] for i in range(len(x) - 1)]
    convex = sum(all(a <= b for a, b in zip(inc(x), inc(x)[1:])) for x in F12)
    concave = sum(all(a >= b for a, b in zip(inc(x), inc(x)[1:])) for x in F12)
    res["k12_convex_concave"] = {"ray_facets": len(F12), "convex": convex, "concave": concave,
                                 "not_convex": len(F12) - convex, "neither": len(F12) - convex - concave}
    checks["k12_convex_order"] = report("k = 12: 7,008 of 7,018 ray facets not convex, 6,998 neither convex nor concave",
                                        len(F12) == 7018 and len(F12) - convex == 7008 and len(F12) - convex - concave == 6998,
                                        f"({t()} s)")

    g = (1, 1, 5, 3, 2, 2, 3, 5); Eg = ST.ends(g); M = Mfun(Eg); att = attain_matrix(Eg, M)
    partners4 = [q for q in range(len(g)) if att[4][q]]
    ok = Eg == [0, 1, 6, 9, 11, 13, 16, 21] and partners4 == [4, 7] and M[8] == 22
    ok &= attainers(Eg, 8) == [(1, 7), (2, 6), (3, 5), (4, 4)]
    cg = classify(ST.u_of(g)); ok &= cg["facet"] and cg["total_support"] and sum(g) == 22
    checks["lemma_6_6_example"] = report("g = (1,1,5,3,2,2,3,5): E, partners of block 4, M(8) = 22 attained by four pairs, "
                                         "total support, facet", ok)

    # ---- total support of staircases
    g = (1, 2, 2, 1, 2, 1, 2, 3, 1, 1, 3, 2, 1, 2, 1, 2, 2, 4, 2, 2); r = len(g) - 1; Eg = ST.ends(g); Sg = ST.starts(g)
    T = tight_edges(ST.u_of(g)); att = attain_matrix(Eg)
    nb7 = sorted(j for (i, j) in T if i == Sg[7])
    A = [[l for l in range(r + 1) if att[q][l]] for q in range(r)]
    Ks = _largest_balanced(A, r)
    Astar = set().union(*[set(A[q]) for q in Ks]); Bstar = {d for d in range(1, r + 1) if g[d] == 1 and d - 1 in Ks}
    rigid = [(s, attainers(Eg, s)[0]) for s in range(2 * r + 1) if len(attainers(Eg, s)) == 1]
    cert = [s for s, (c0, d0) in rigid if (c0 in Astar and d0 not in Bstar) or (d0 in Astar and c0 not in Bstar)]
    ok = sum(g) == 37 and Sg[7] == 11 and Eg[10] == 18 and nb7 == [18] and A[6] == [10]
    ok &= not any(10 in pair for _, pair in rigid)
    ok &= sorted(q + 1 for q in Ks) == [6, 7, 10] and len(cert) == 9
    res["minimal_balanced_example"] = {"neighbors_of_S7": nb7, "K_star_starts": sorted(q + 1 for q in Ks),
                                       "certifying_rigid_levels": cert}
    checks["minimal_balanced_example"] = report("k = 37: {S_7} minimal balanced with N = {E_10}, no rigid level through 10, "
                                                "K* = {S_6, S_7, S_10} certified by nine rigid levels", ok)

    g = (1, 1, 4, 4, 7, 4, 4, 1, 6, 5); Eg = ST.ends(g); u = ST.u_of(g)
    cl = boundary_classes(Eg); sol = solitary_indices(Eg); att = attain_matrix(Eg)
    A7 = [l for l in range(len(g)) if att[7][l]]
    ok = sum(g) == 37 and cl == [[0, 1, 2, 3, 4, 5, 6, 7], [8]] and sol == [7] and A7 == [0]
    ok &= perfect_matching(len(u), tight_edges(u)) is not None and not has_total_support(u)
    res["solitary_example"] = {"classes": cl, "solitary": sol, "A(7)": A7}
    checks["solitary_example"] = report("k = 37: classes {E_0..E_7}, {E_8}; only solitary index 7, A(7) = {0}; "
                                        "tight permutation, no total support", ok)
    # ---- small k: coordinate facets and facet rays (after tab:counts)
    from dsxray import small_k_fan as SK
    ok = True; small = {}
    for k in (2, 3, 4):
        H, _ = SK.exact_facets(k)
        coord = {SK.coordinate_inequality(k, tt) for tt in range(2 * k - 1) if classify_coordinate(k, tt)["facet"]}
        steps = {SK.ray_inequality([s * int(i >= j) for i in range(k)]) for j in range(1, k) for s in (1, -1)}
        small[k] = {"facets": len(H), "coordinate_facets": len(coord & H), "ray_facets": len(H - coord)}
        ok &= (H - coord) <= steps and coord <= H
    ok &= [small[k]["facets"] for k in (2, 3, 4)] == [2, 6, 13]
    ok &= [small[k]["coordinate_facets"] for k in (2, 3, 4)] == [2, 4, 7]
    res["small_k"] = small
    checks["small_k_facets"] = report("k = 2, 3, 4: 2, 4, 7 coordinate facets; every other facet is a signed single step",
                                      ok, str(small))

    # ---- app:computation: deterministic first facet of the walk
    def face_dim(u):
        kk = len(u); rows = []
        for tt in range(2 * kk - 1):
            P = [(i, tt - i) for i in range(kk) if i <= tt - i < kk]
            m = min(u[i] + u[j] for i, j in P); Smin = [(i, j) for i, j in P if u[i] + u[j] == m]
            def vec(pr):
                v = [0] * kk; v[pr[0]] += 1; v[pr[1]] += 1
                return v
            b0 = vec(Smin[0])
            rows += [[x - y for x, y in zip(vec(pr), b0)] for pr in Smin[1:]]
        return rank_exact(rows) if rows else 0
    ok = all(face_dim([max(0, i - 1) for i in range(k)]) == k - 3 for k in range(4, 16))
    checks["hinge_first_facet"] = report("app:computation: dim Q^{h_1} = k - 3 for 4 <= k <= 15 (exact rank)", ok)

    # ---- sec:s_census-proofs: a ray with two classes of boundary steps
    g = (1, 2, 1, 3); cl = boundary_classes(ST.ends(g))
    ok = classify(ST.u_of(g))["ray"] and len(cl) == 2
    checks["ray_two_classes"] = report("g = (1, 2, 1, 3) (k = 7) is a ray with two classes of boundary steps", ok, str(cl))

    res["checks"] = checks
    finish(res, out_dir(a.out), "examples")


if __name__ == "__main__":
    main()
