#!/usr/bin/env python3
"""ex:noTS and rem:samefacet (a facet ray without total support at k = 19; two rays, one facet), and further facet
rays and potentials without total support (data/counterexamples/potentials.json).

  1. ex:noTS, chi = indicator of O = {2,4,5,7,9,10,11,12,13} on {0..18}: the sums of two elements of O^c cover
     {0..36} except 5, 10, 13, 27, which are the stated sums of an element of O^c and an element of O; v = 1 exactly
     on {5, 10, 13, 27}; the only tight neighbor of 11 is 16; the tight pairs outside every tight permutation are
     exactly the 19 tight pairs that contain 16 and not 11; chi is a ray (dim Z = 3, exact); every anti-diagonal
     carries an allowed edge; the face has dimension 34 = 2k - 4 (exact rank, and the certified facet test); {11, 16}
     is a bipartite component of the allowed graph, bip = 1, dim Z* = 4; chi + iota has increments in {0, 1, 2},
     height two.
  2. rem:samefacet, chi' = chi - e_11 + e_16: the same v and the same sum as chi; a ray; the only tight neighbor of 16
     is 11; chi' is not in span{1, iota, chi} (exact rank 4); chi - chi' lies in Y*(chi); I(chi') defines a facet.
  3. Two potentials at k = 17 (not cited in the paper): W, not a ray, whose tight graph has total support and the
     bipartite component {8, 10}, with W + eps (e_8 - e_10) having the same v and the same tight graph for
     -1 < eps < 1 (checked at eps = +-1/2, +-9/10, new tight pairs at eps = +-1); and the halved potential at
     eps = 1, which has a tight permutation, no total support, an allowed edge on every anti-diagonal, and gives no
     facet.
  4. The family R(m) = (0,0,2,2,0,0,1,1,1,2,2,0,0,2,2,0,0, 2^m) (k = 17 + m; not cited in the paper): R(3), R(4) are
     not rays, and R(m) is a facet ray without total support for 5 <= m <= --rmax (default 7; --full: 24).
  5. With --search: rem:samefacet, no facet ray of height at most two without total support for 5 <= k <= --kmax
     (alias --search-kmax; default 12, seconds; --full sets 18, the paper's range; on one core k <= 16 takes about
     5 minutes, k = 17 about 6 minutes and k = 18 about 18 minutes, about 28 minutes in all), by the C program of
     scripts/search_total_support.py with ranks modulo two primes (every reported potential re-verified exactly).
Default: items 1-4 (about 10 s).
"""
import argparse
import os
from fractions import Fraction

from common import DATA, Timer, build_c, finish, out_dir, read_json, report

from dsxray.facet_criterion import allowed_graph, components, criterion, rank_int
from dsxray.facet_test import classify, tight_edges, v_of
from dsxray.padding import pattern

from search_total_support import run_increments


def potential(name):
    for p in read_json(os.path.join(DATA, "counterexamples", "potentials.json"))["potentials"]:
        if p["name"] == name:
            return p["u"]
    raise KeyError(name)


def R(m):
    return [0, 0, 2, 2, 0, 0, 1, 1, 1, 2, 2, 0, 0, 2, 2, 0, 0] + [2] * m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/facets_without_total_support")
    ap.add_argument("--rmax", type=int, default=7)
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--search-kmax", "--kmax", dest="search_kmax", type=int, default=12)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 24); a.search = True; a.search_kmax = max(a.search_kmax, 18)
    out = out_dir(a.out); res = {}; checks = {}

    # 1. ex:noTS
    t = Timer(); k = 19; O = {2, 4, 5, 7, 9, 10, 11, 12, 13}; Oc = sorted(set(range(k)) - O)
    chi = potential("example_5_1_chi"); ok = chi == [1 if i in O else 0 for i in range(k)]
    missing = sorted(set(range(2 * k - 1)) - {a_ + b_ for a_ in Oc for b_ in Oc})
    ok &= missing == [5, 10, 13, 27]
    ok &= all(x in O and y in Oc and x + y == s_ for (x, y), s_ in zip([(5, 0), (10, 0), (13, 0), (11, 16)], missing))
    v = v_of(chi); ok &= [t_ for t_ in range(2 * k - 1) if v[t_] == 1] == missing and set(v) == {0, 1}
    E, A, _ = allowed_graph(chi)
    ok &= sorted(j for j in range(k) if (11, j) in E) == [16]
    lost = E - A
    ok &= lost == {(i, j) for (i, j) in E if 16 in (i, j) and 11 not in (i, j)} and len(lost) == 19
    c = criterion(chi, face_dim=True)
    ok &= c["zdim"] == 3 and c["cover"] and c["face_dim"] == 34 == 2 * k - 4 and classify(chi)["facet"]
    ok &= c["bip"] == 1 and ([11, 16], True) in components(k, A) and c["zdim_star"] == 4
    p = pattern([chi[i] + i for i in range(k)])
    inc = [chi[i + 1] + 1 - chi[i] for i in range(k - 1)]
    ok &= set(inc) == {0, 1, 2} and max(p) == 2 and p == tuple(inc)
    res["example_5_1"] = {k_: c[k_] for k_ in ("zdim", "lost_edges", "cover", "face_dim", "bip", "zdim_star")}
    checks["example_5_1"] = report("ex:noTS: sums, v, the unique neighbor 16 of 11, 19 lost tight pairs, ray, facet (dimension 34), bip = 1, dim Z* = 4, height two", ok, f"({t()} s)")

    # 2. rem:samefacet
    t = Timer(); chp = potential("remark_5_2_chi_prime")
    ok = chp == [x - (i == 11) + (i == 16) for i, x in enumerate(chi)]
    ok &= v_of(chp) == v and sum(chp) == sum(chi)
    cp = criterion(chp); ok &= cp["zdim"] == 3 and classify(chp)["facet"]
    Ep = tight_edges(chp); ok &= sorted(j for j in range(k) if (16, j) in Ep) == [11]
    ok &= rank_int([[1] * k, list(range(k)), chi, chp]) == 4
    d = [x - y for x, y in zip(chi, chp)]; ok &= all(d[i] + d[j] == 0 for (i, j) in A)
    res["remark_5_2"] = {"lost_edges_chi_prime": cp["lost_edges"]}
    checks["remark_5_2"] = report("rem:samefacet: chi' has the same v and sum, is a ray giving a facet, unique neighbor 11 of 16, not in span{1, iota, chi}, chi - chi' in Y*(chi)", ok, f"({t()} s)")

    # 3. k = 17
    t = Timer(); W = potential("k17_total_support_bipartite_component"); h = potential("k17_cover_without_total_support")
    k = 17; cW = criterion(W); EW, AW, _ = allowed_graph(W)
    ok = cW["total_support"] and cW["zdim"] == 6 and ([8, 10], True) in components(k, EW)
    vW = v_of(W)
    for eps in (Fraction(1, 2), Fraction(-1, 2), Fraction(9, 10), Fraction(-9, 10)):
        We = [W[i] + eps * ((i == 8) - (i == 10)) for i in range(k)]
        ok &= v_of(We) == vW and tight_edges(We) == EW
    for eps in (1, -1):
        We = [W[i] + eps * ((i == 8) - (i == 10)) for i in range(k)]
        ok &= v_of(We) == vW and tight_edges(We) > EW
    ok &= h == [x // 2 for x in [W[i] + (i == 8) - (i == 10) for i in range(k)]]
    ch = criterion(h, face_dim=True)
    ok &= ch["tight_permutation"] and not ch["total_support"] and ch["cover"] and not classify(h)["facet"]
    res["k17"] = {"W": {x: cW[x] for x in ("zdim", "total_support", "bip")}, "halved": {x: ch[x] for x in ("zdim", "lost_edges", "cover", "face_dim")}}
    checks["k17"] = report("k = 17: total support with a bipartite component {8, 10} (not a ray); the halved potential has an allowed edge on every anti-diagonal without total support", ok, f"({t()} s)")

    # 4. R(m)
    t = Timer(); ok = True; fam = {}
    for m in range(3, a.rmax + 1):
        u = R(m); c = criterion(u); f = classify(u)["facet"]
        fam[m] = {"k": len(u), "zdim": c["zdim"], "total_support": c["total_support"], "facet": f}
        ok &= (c["zdim"] > 3 and not f) if m < 5 else (c["zdim"] == 3 and not c["total_support"] and f and c["criterion_facet"])
    ok &= R(5) == potential("k22_R5")
    res["R_family"] = fam
    checks["R_family"] = report(f"R(m): not rays for m = 3, 4; facet rays without total support for 5 <= m <= {a.rmax}", ok, f"({t()} s)")

    # 5. search
    if a.search:
        t = Timer(); exe = build_c("search_total_support", os.path.join(out, "bin")); found = []; tab = {}
        for kk in range(5, a.search_kmax + 1):
            tot, f = run_increments(exe, kk, 2, a.workers); tab[kk] = tot; found += f
            print(f"   k={kk}: {tot['rays_with_tight_permutation']} rays of height <= 2 with a tight permutation, {tot['rays_without_total_support']} without total support, {len(f)} facets ({t()} s)", flush=True)
        res["search_height_2"] = {"table": tab, "found": found}
        checks["search_height_2"] = report(f"rem:samefacet: no facet ray of height <= 2 without total support for k <= {a.search_kmax}",
                                           not found, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "facets_without_total_support")


if __name__ == "__main__":
    main()
