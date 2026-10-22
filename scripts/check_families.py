#!/usr/bin/env python3
"""Infinite families of facets from no step potential (thm:fam and sec:SM-famproof), exact.

Notation of thm:fam: K = k-1, L = k-2; (a) u_{a,b}(i) = min(i, a) on [0, n], u_{a,b}(n+1) = u_{a,b}(n+2) = a+2
(a >= 2, b >= 1, n = a+b, k = n+3); (b) omega = e_1 - e_L; (c) nu = e_K - e_0 - e_2. In the code: family_A_u,
W_u (omega), Vplus_u (nu), Vminus_u (-nu).
  1. Membership (certified facet test, exact ray test, total support): u_{a,b} for 6 <= k <= --a-kmax; omega and
     -omega for k = 7 and 9 <= k <= --kmax; nu for 7 <= k <= --kmax; -nu for k = 7 and 9 <= k <= --kmax.
     sec:SM-famproof: omega and -nu are not rays for k = 6 and k = 8, nu is not a ray for k = 6, and
     -u_{a,b} defines no facet for 6 <= k <= --a-kmax (the paper states 6 <= k <= 12).
  2. The formulas for v = u box u in the proof, and the tightness of every pair the proof uses (the attaining pairs,
     the universal index e, the tight permutation fixing e, and the pairs of the paragraphs computing Z(u)), for
     7 <= k <= --proof-kmax (omega, nu, -nu) and 6 <= k <= --a-kmax (u_{a,b}); the pairs of the Z(u) paragraphs are
     checked where those paragraphs apply (k >= 9 for omega and -nu, k >= 7 for nu). Also the complete tight graphs
     predicted in the code agree with the computed ones.
  3. thm:fam(e): P_L and P_R hold for every member with 9 <= k <= --proof-kmax, so all paddings are facets.
  4. Coverage: the families, their reflections and paddings give 44, 80, 130, 196 of the 195, 622, 1,902, 5,748
     non-step facets for k = 9, 10, 11, 12 (complete facet lists, data/facets).
Defaults: --kmax 16, --a-kmax 16, --proof-kmax 24. --full: 30, 22, 30.
"""
import argparse

from common import Timer, finish, load_facet_rows, out_dir, report

from dsxray import nonstep_families as NF
from dsxray.facet_test import classify, is_ray, tight_edges, v_of
from dsxray.padding import P_L, P_R, is_step, pattern

COVER = {9: (195, 44), 10: (622, 80), 11: (1902, 130), 12: (5748, 196)}


def proof_data(name, k, a=None):
    """(v predicted by the proof, list of pairs the proof uses where its paragraphs apply)."""
    K, L = k - 1, k - 2; T = 2 * k - 1; pairs = []
    if name == "u":
        b = k - 3 - a; n = a + b; v = []
        for t in range(T):
            if t <= n:
                v.append(min(t, a))
            elif t <= 2 * n:
                v.append(a + min(t - n - 1 + 1, a))
            else:
                v.append(2 * a + 2 if t <= 2 * n + 2 else 2 * a + 4)
        pairs += [(n, j) for j in range(k)]
        pairs += [(s - 1, n + 2) for s in range(1, n) if s + 1 <= a]
        pairs += [(n, n + 1), (n - 1, n + 2), (n, n + 2)]
        pairs += [p for j in range(1, a) for p in ((0, j + 1), (1, j))]
        pairs += [p for i in range(a, n) for p in ((i, n), (i + 1, n - 1))]
        pairs += [(a, n), (a - 2, n + 2)]
        perm = [a - i for i in range(a + 1)] + [a + n - i for i in range(a + 1, n)] + [n, n + 2, n + 1]
        e = n
    elif name == "omega":
        v = [0, 1] + [0] * (L - 2) + [-1, 0] + [-1] * (L - 2) + [-2, -1, 0]
        Q = range(2, L)
        pairs += [(0, t) for t in range(2, L)] + [(i, j) for i in Q for j in Q if i + j <= L - 1]
        pairs += [(0, L), (0, K), (1, L)] + [(i, L + 1 - i) for i in Q] + [(t - L, L) for t in range(L + 2, 2 * L)]
        if k >= 9:
            pairs += [p for t in range(6, L) for p in ((2, t - 2), (3, t - 3))]
            pairs += [(2, L - 1), (3, L - 2), (4, L - 3), (0, 4), (2, 2)]
        perm = [1, 0] + [L + 1 - i for i in Q] + [L, K]; e = L
    elif name == "nu":
        v = [-2, -1, -2, -1, -2] + [-1] * (L - 2) + [0] * (L - 2) + [1, 2]
        pairs += [(0, 0), (0, 2), (2, 2), (2, 1)] + [(2, t - 2) for t in range(5, L + 3)] + [(2, K)]
        pairs += [(i, j) for i in range(3, L + 1) for j in range(3, L + 1) if i + j >= L + 3]
        if k >= 7:
            pairs += [p for i in range(3, L) for p in ((i, L), (i + 1, L - 1))]
            pairs += [(0, 5), (2, 3), (0, 3), (1, 2), (3, L)]
        perm = [1, 0, 2] + [L + 3 - i for i in range(3, L + 1)] + [K]; e = 2
    elif name == "-nu":
        v = [2, 1, 0, 1] + [0] * (L - 2) + [-1, 0] + [-1] * (L - 2) + [-2]
        Z0 = [1] + list(range(3, L + 1))
        pairs += [(i, j) for i in Z0 for j in Z0 if 4 <= i + j <= L + 1] + [(0, K), (1, K), (2, K)]
        pairs += [(t - K, K) for t in range(L + 4, 2 * L + 2)]
        if k >= 9:
            pairs += [p for t in range(8, L + 2) for p in ((3, t - 3), (4, t - 4))]
            pairs += [(3, L), (4, L - 1), (5, L - 2), (1, 5), (3, 3), (0, 3), (1, 2), (1, L)]
        perm = [0, 2, 1] + [L + 3 - i for i in range(3, L + 1)] + [K]; e = K
    pairs += [(e, j) for j in range(k)] + [(i, perm[i]) for i in range(k)]
    return v, pairs, perm


def potential(name, k, a=None):
    if name == "u":
        return NF.family_A_u(a, k - 3 - a)
    return {"omega": NF.W_u(k), "nu": NF.Vplus_u(k), "-nu": NF.Vminus_u(k)}[name]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/families")
    ap.add_argument("--kmax", type=int, default=16)
    ap.add_argument("--a-kmax", type=int, default=16)
    ap.add_argument("--proof-kmax", type=int, default=24)
    ap.add_argument("--cover-kmax", type=int, default=12)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 30); a.a_kmax = max(a.a_kmax, 22); a.proof_kmax = max(a.proof_kmax, 30)
    res = {}; checks = {}

    t = Timer(); bad = []
    for k in range(6, a.a_kmax + 1):
        for aa in range(2, k - 3):
            u = NF.family_A_u(aa, k - 3 - aa)
            if NF.check_member(u) != (True, True, True):
                bad.append(("u_ab", k, aa))
            if classify([-x for x in u])["facet"]:
                bad.append(("-u_ab is a facet", k, aa))
    for k in range(6, a.kmax + 1):
        exp_w = k == 7 or k >= 9; exp_vm = k == 7 or k >= 9; exp_vp = k >= 7
        for name, u, exp in (("omega", NF.W_u(k), exp_w), ("-omega", [-x for x in NF.W_u(k)], exp_w),
                             ("nu", NF.Vplus_u(k), exp_vp), ("-nu", NF.Vminus_u(k), exp_vm)):
            f, ray, ts = NF.check_member(u)
            if f != exp or (exp and not (ray and ts)):
                bad.append((name, k))
    res["membership_exceptions"] = bad
    checks["membership"] = report(f"u_ab (k <= {a.a_kmax}), omega, -omega, nu, -nu (k <= {a.kmax}) are facets "
                                  f"exactly where stated; -u_ab never", not bad, f"({t()} s)")
    nonray = {"omega": [is_ray(NF.W_u(k)) for k in (6, 8)], "-nu": [is_ray(NF.Vminus_u(k)) for k in (6, 8)],
              "nu": [is_ray(NF.Vplus_u(6))]}
    res["not_rays"] = nonray
    checks["not_rays"] = report("omega and -nu are not rays for k = 6, 8; nu is not a ray for k = 6",
                                not any(x for v in nonray.values() for x in v))

    t = Timer(); vbad = []; pbad = []
    jobs = [("u", k, aa) for k in range(6, a.a_kmax + 1) for aa in range(2, k - 3)]
    jobs += [(nm, k, None) for k in range(7, a.proof_kmax + 1) for nm in ("omega", "nu", "-nu")]
    for nm, k, aa in jobs:
        u = potential(nm, k, aa); v, pairs, _ = proof_data(nm, k, aa); E = tight_edges(u)
        if v_of(u) != v:
            vbad.append((nm, k, aa))
        if not all(0 <= i < k and 0 <= j < k and (i, j) in E for i, j in pairs):
            pbad.append((nm, k, aa))
    res["v_formula_exceptions"] = vbad; res["pair_exceptions"] = pbad
    checks["v_formulas"] = report(f"formulas for v in the proof (omega, nu, -nu: 7 <= k <= {a.proof_kmax}; "
                                  f"u_ab: 6 <= k <= {a.a_kmax})", not vbad, f"({t()} s)")
    checks["proof_pairs"] = report("every pair used in the proof is tight", not pbad)
    ok = True
    for k in range(9, a.proof_kmax + 1):
        ok &= tight_edges(NF.W_u(k)) == NF.W_tight_graph(k) and NF.permutation_with_universal_vertex(NF.W_tight_graph(k), *NF.W_permutation(k))
        ok &= tight_edges(NF.Vplus_u(k)) == NF.Vplus_tight_graph(k) and NF.permutation_with_universal_vertex(NF.Vplus_tight_graph(k), *NF.Vplus_permutation(k))
        ok &= tight_edges(NF.Vminus_u(k)) == NF.Vminus_tight_graph(k) and NF.permutation_with_universal_vertex(NF.Vminus_tight_graph(k), *NF.Vminus_permutation(k))
    for k in range(6, a.a_kmax + 1):
        for aa in range(2, k - 3):
            b = k - 3 - aa; E = tight_edges(NF.family_A_u(aa, b))
            ok &= E == NF.family_A_tight_graph(aa, b) and NF.permutation_with_universal_vertex(E, *NF.family_A_permutation(aa, b))
    checks["tight_graphs"] = report("complete tight graphs, universal index and tight permutation as predicted", ok)

    t = Timer(); ok = True
    for k in range(9, a.proof_kmax + 1):
        us = [NF.W_u(k), [-x for x in NF.W_u(k)], NF.Vplus_u(k), NF.Vminus_u(k)]
        us += [NF.family_A_u(aa, k - 3 - aa) for aa in range(2, k - 3)]
        ok &= all(P_L(pattern(u)) and P_R(pattern(u)) for u in us)
    checks["padding_conditions"] = report(f"(e): P_L and P_R for all members, 9 <= k <= {a.proof_kmax}", ok, f"({t()} s)")

    t = Timer(); cover = {}; ok = True
    for k in range(9, a.cover_kmax + 1):
        F = {pattern(r["u"]) for r in load_facet_rows(k) if r["facet"]}
        NS = {p for p in F if not is_step(p)}
        P = NF.padding_closure(NF.proved_patterns, k)
        ok &= P <= F
        row = (len(NS), len(P & NS)); cover[k] = dict(zip(("nonstep_facets", "from_theorem_5_3"), row))
        ok &= row == COVER[k]
    res["coverage"] = cover
    checks["coverage"] = report("coverage of the non-step facets, k = 9..12: 44, 80, 130, 196", ok, f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "families")


if __name__ == "__main__":
    main()
