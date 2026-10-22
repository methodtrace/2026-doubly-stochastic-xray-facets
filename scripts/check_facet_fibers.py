#!/usr/bin/env python3
"""Facet fibers: the fixed-coefficient fiber of a facet, and the cubical fibers on 19m points.

Claims (thm:fiber and prop:cubical, proved in the paper):
  (a) thm:fiber: for a facet F = I(u_0) (not a coordinate facet) with v_0 = u_0 box u_0, A = E*(u_0) and
      b = bip(A), the fiber D(F) = {w : w(i) + w(j) >= v_0(i + j), 2 sum w = 2 sum u_0} is a compact polytope of
      dimension b with affine hull u_0 + Y*(u_0); every w in D(F) has allowed graph A; w is in the relative interior
      iff E(w) has total support; dim Z(w) = 3 + bip(E(w)); the vertices are exactly the normalized rays defining F.
      So b = 0 gives exactly one ray, with total support; b >= 1 gives at least b + 1 rays, none with total support;
  (b) prop:cubical (with f = m): for chi_m(i) = chi(floor(i/m)) on k = 19m points (chi the potential of ex:noTS),
      v = 1 exactly on {6m - 1, 11m - 1, 14m - 1, 28m - 1}; 19 m^2 tight edges are not allowed; dim Z(A) = 3 + m,
      bip(A) = m, face dimension 2k - 4: a facet without total support; its fiber is the cube
      chi_m + sum_a lambda_a (e_{11m+a} - e_{17m-1-a}), lambda in [-1, 0]^m, with 2^m vertices (2^m distinct rays
      defining one facet; the potential cone is not simplicial for m >= 2).
Checks (exact; allowed edges by forcing each tight edge into a perfect matching, ranks over Q):
  1. (b) for m = 1..--mmax: all statements above, the 2^m vertices are rays with the same v, the centre of the cube
     has tight graph exactly A; for m <= 2 also with facet_criterion.criterion (switching digraph, other ranks);
  2. (a) on the cube: at every point lambda in {-1, -1/2, 0}^m, dim Z(E(w)) = 3 + bip(E(w)), E(w) contains A, and
     E(w) = A (total support) exactly at the centre; the zero set S of chi has S + S = [0, 36] minus {5, 10, 13, 27};
  3. (a) for b = 0: distinct facet rays with total support in data/facets (8 <= k <= --kmax, normalized
     u(0) = u(k-1) = 0) have distinct v, so each defines its own facet.
Default: m <= 2, k <= 10 (about 5 s). --full: m <= 3, k <= 12 (about 30 s).
"""
import argparse
from fractions import Fraction
from itertools import product

from common import Timer, finish, load_facet_rows, out_dir, report

from dsxray import exact_graphs as GR
from dsxray import facet_criterion as FC
from dsxray import potential_families as PF


def cube_points(m):
    u = PF.chi(m); k = len(u)
    _, tight = GR.tight_graph(u); A = GR.allowed_edges(k, tight)
    pairs = [(11 * m + a, 17 * m - 1 - a) for a in range(m)]
    ok = True
    for lam in product([Fraction(-1), Fraction(-1, 2), Fraction(0)], repeat=m):
        w = [Fraction(x) for x in u]
        for (i, j), x in zip(pairs, lam):
            w[i] += x; w[j] -= x
        _, E = GR.tight_graph(w)
        bip = sum(c["bipartite"] for c in GR.components(k, E))
        ok &= GR.zdim(k, E) == 3 + bip and A <= E
        ok &= (E == A) == all(x == Fraction(-1, 2) for x in lam)
    return ok, 3 ** m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/facet_fibers")
    ap.add_argument("--mmax", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.mmax = max(a.mmax, 3); a.kmax = max(a.kmax, 12)
    out = out_dir(a.out); res = {}; checks = {}
    for m in range(1, a.mmax + 1):
        t = Timer()
        try:
            s = PF.cubical_fiber(m); ok = True
        except AssertionError as e:
            s, ok = {"error": str(e)}, False
        if ok and m <= 2:
            cr = FC.criterion(PF.chi(m), face_dim=True)
            ok &= (cr["zdim"], cr["zdim_star"], cr["bip"], cr["face_dim"], cr["lost_edges"], cr["total_support"]) == \
                (3, 3 + m, m, 38 * m - 4, 19 * m * m, False)
        res[f"cube_m{m}"] = s
        checks[f"cube_m{m}"] = report(f"cubical fiber m = {m}: facet without total support, 2^{m} rays defining it", ok,
                                      f"{s} ({t()} s)")
        t = Timer(); ok, n = cube_points(m)
        checks[f"fiber_points_m{m}"] = report(f"m = {m}: dim Z(w) = 3 + bip(E(w)) on the cube, total support only at the centre",
                                              ok, f"{n} points ({t()} s)")
    checks["zero_set"] = report("S + S = [0, 36] minus {5, 10, 13, 27} for the zero set S of chi", PF.zero_set_certificate())
    t = Timer(); ok = True; n = 0
    for k in range(8, a.kmax + 1):
        seen = {}
        for row in load_facet_rows(k):
            if row["facet"] and row["total_support"]:
                v = tuple(GR.tight_graph(row["u"])[0]); n += 1
                ok &= v not in seen
                seen[v] = row["u"]
    res["unique_rays"] = {"k_range": [8, a.kmax], "facet_rays_with_total_support": n}
    checks["unique_rays"] = report(f"b = 0: distinct facet rays with total support have distinct v, 8 <= k <= {a.kmax}", ok,
                                   f"{n} rays ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "facet_fibers")


if __name__ == "__main__":
    main()
