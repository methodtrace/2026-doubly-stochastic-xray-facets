#!/usr/bin/env python3
"""Sparse refinement in reverse (descent), its counterexample without total support, and ray-only padding.

Claims ((a)-(c) are prop:descent; (d) is thm:pad(b) together with the ray part of the proof of thm:pad(c), as used in
sec:SM-padcount), for a normalized potential u on {0..K}, m >= 2 and u' = spr_m u (u'(i) = u(ceil(i/m))):
  (a) if E(u') has total support, so does E(u); if the bipartite tight graph of u' is connected, so is that of u
      (no ray hypothesis): for a base tight edge (a, b) take a refined tight permutation through (ma, mb); its block
      counts C (blocks {0}, {m(a-1)+1..ma}) give the doubly stochastic D = C/m + (1 - 1/m) e_0 e_0^T supported on E(u)
      with D_ab > 0;
  (b) if moreover u' is a ray and the flat-level common-difference condition holds (the values h_z(s) - h_z(s-1),
      s in F, agree for every z in Z(u); automatic when u has exactly one zero increment), then u is a ray;
  (c) without total support the ray conclusion fails: for increments (0,1,1,2,0) and m = 2 (also m = 3), spr_m u is
      a ray and u is not (dim Z(u) = 4); the condition of (b) fails there and spr_m u lacks total support;
  (d) ray-only padding: for a normalized ray u, pad_L u is a ray iff P_L holds (an index q with d(q) = 0 and (0, q)
      tight), and pad_R u is a ray iff P_R holds; no total support or connectivity needed; P_L persists under
      repeated left padding.
Checks (exact; matchings by forcing edges, ranks over Q):
  1. (a), (b) every normalized word with increments in {0, 1, 2}, k <= --kmax points, m = 2, 3: projection and lift
     of tight edges, descent of total support and connectivity, every collapse certificate (refined permutation,
     row and column sums m, support in E(u), positive at the edge), and the ray descent under the condition;
     total support also compared with facet_test.has_total_support;
  2. (c) the six words (0,1,1,2,0) (m = 2, 3), (0,0,1,1,2,0), (0,1,1,1,2,0), (0,1,1,2,0,2), (1,0,1,1,2,0) (m = 2);
     dim Z by facet_test.z_dim as well;
  3. (d) every primitive normalized word with increments in {0, 1, 2} and k <= --pad-kmax points: the two
     equivalences for every ray, and pad_L^j (j = 1, 2, 3) of rays with P_L are rays with P_L.
--kmax 5: 172 base/factor pairs; --pad-kmax 7: 846 words, 171 rays.
Default: k <= 5, padding k <= 6 (about 2 s). --full: k <= 6, padding k <= 7 (about 5 s).
"""
import argparse
from functools import reduce
from itertools import product
from math import gcd

from common import Timer, finish, out_dir, report

from dsxray import exact_graphs as GR
from dsxray import potential_families as PF
from dsxray.facet_test import has_total_support, z_dim

COUNTEREXAMPLES = [((0, 1, 1, 2, 0), 2), ((0, 1, 1, 2, 0), 3), ((0, 0, 1, 1, 2, 0), 2), ((0, 1, 1, 1, 2, 0), 2),
                   ((0, 1, 1, 2, 0, 2), 2), ((1, 0, 1, 1, 2, 0), 2)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/refinement_descent")
    ap.add_argument("--kmax", type=int, default=5)
    ap.add_argument("--pad-kmax", type=int, default=6)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 6); a.pad_kmax = max(a.pad_kmax, 7)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); ok = True; c = {"pairs": 0, "refined_total_support": 0, "certificates": 0, "refined_rays": 0,
                                 "condition_holds_on_refined_rays": 0}
    for k in range(3, a.kmax + 1):
        for d in product(range(3), repeat=k - 1):
            if min(d) != 0 or max(d) == 0:
                continue
            u = GR.from_increments(d)
            for m in (2, 3):
                try:
                    bts, rts, bc, rc, bz, rz, cond, nc = PF.descent(u, m)
                except AssertionError:
                    ok = False; continue
                ok &= bts == has_total_support(u) and rts == has_total_support(PF.spr(u, m))
                c["pairs"] += 1; c["refined_total_support"] += rts; c["certificates"] += nc
                c["refined_rays"] += rz == 3; c["condition_holds_on_refined_rays"] += rz == 3 and cond
    res["descent"] = dict(c, kmax=a.kmax)
    if a.kmax == 5:
        ok &= c["pairs"] == 172
    checks["descent"] = report(f"descent of total support and connectivity, ray descent under the condition, k <= {a.kmax}, m = 2, 3",
                               ok, f"{c} ({t()} s)")

    t = Timer(); ok = True; rows = []
    for d, m in COUNTEREXAMPLES:
        u = GR.from_increments(d); up = PF.spr(u, m)
        v, E = GR.tight_graph(u); _, Ep = GR.tight_graph(up)
        bz, rz = GR.zdim(len(u), E), GR.zdim(len(up), Ep)
        good = bz == 4 == z_dim(u) and rz == 3 == z_dim(up) and not PF.lifting_condition(u, v, E)
        good &= not PF.total_support(len(up), Ep)
        ok &= good; rows.append({"increments": list(d), "m": m, "base_zdim": bz, "refined_zdim": rz})
    res["counterexamples"] = rows
    checks["counterexamples"] = report("without total support ray descent fails: (0,1,1,2,0), m = 2, 3, and four more words",
                                       ok, f"({t()} s)")

    t = Timer(); ok = True; c = {"words": 0, "rays": 0, "P_L": 0, "P_R": 0, "iterated_left": 0}
    for k in range(3, a.pad_kmax + 1):
        for w in product(range(3), repeat=k - 1):
            if min(w) != 0 or max(w) == 0 or reduce(gcd, w) != 1:
                continue
            c["words"] += 1
            try:
                ray, pl, pr = PF.padding_ray_only(w)
            except AssertionError:
                ok = False; continue
            if not ray:
                continue
            c["rays"] += 1; c["P_L"] += pl; c["P_R"] += pr
            if pl:
                u = GR.from_increments(w)
                for _ in range(3):
                    u = [u[0]] + u
                    _, E = GR.tight_graph(u)
                    ok &= GR.zdim(len(u), E) == 3 and 0 in PF.padding_witnesses(u, E)[0]
                    c["iterated_left"] += 1
    res["padding"] = dict(c, kmax=a.pad_kmax)
    if a.pad_kmax == 7:
        ok &= (c["words"], c["rays"]) == (846, 171)
    checks["padding"] = report(f"ray-only padding: pad_L ray iff P_L, pad_R ray iff P_R, k <= {a.pad_kmax}", ok, f"{c} ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "refinement_descent")


if __name__ == "__main__":
    main()
