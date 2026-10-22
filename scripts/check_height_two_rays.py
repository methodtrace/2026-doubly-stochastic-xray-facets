#!/usr/bin/env python3
"""Height-two rays of the fan: an exponential lower bound from padded cubical corners.

Claim (rem:padcount, proved in sec:SM-padcount): the 2^m vertices of the cubical fiber of chi_m (check_facet_fibers.py),
normalized, are words with increments in {0, 1, 2}, first and last increment 1, with P_L witness q = 14m - 1 and
P_R witness q = 9m - 1; by ray-only padding every pad_L^a pad_R^b of them is a ray, and the words are distinct
(a is the leading run of 0s, b the trailing run of 2s). Hence for k >= 19 the number of height-two rays of the fan
of hat-P(k) (modulo affine functions and positive scaling) is at least
    sum_{m=1}^{floor(k/19)} 2^m (k - 19m + 1)  >=  2^{floor(k/19)}.
(This counts rays, not distinct facets: the 2^m corners define the same facet.)
Checks (exact): for 19 <= k <= --kmax, all padded corner words with 19m + a + b = k are pairwise distinct, normalized
of height 2, rays (dim Z = 3 by rational elimination, and by facet_test.z_dim for k <= 25), and their number equals
the sum above; for m <= 2 the corners have first and last increment 1 and the stated witnesses.
Default: k <= 24 (about 5 s). --full: k <= 40 (about 30 s).
"""
import argparse

from common import Timer, finish, out_dir, report

from dsxray import exact_graphs as GR
from dsxray import potential_families as PF
from dsxray.facet_test import z_dim


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/height_two_rays")
    ap.add_argument("--kmax", type=int, default=24)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 40)
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer(); ok = True
    for m in (1, 2):
        for d in PF.cubical_corners(m):
            u = GR.from_increments(d); _, E = GR.tight_graph(u)
            pl, pr = PF.padding_witnesses(u, E)
            ok &= d[0] == 1 and d[-1] == 1 and max(d) == 2 and min(d) == 0 and 14 * m - 1 in pl and 9 * m - 1 in pr
    checks["corners"] = report("normalized corners: increments in {0,1,2}, first and last 1, witnesses 14m-1 and 9m-1", ok,
                               f"({t()} s)")
    t = Timer(); ok = True; rows = []
    for k in range(19, a.kmax + 1):
        words = PF.padded_corner_words(k)
        ws = [w for _, _, _, w in words]
        good = len(set(ws)) == len(ws) == sum(2 ** m * (k - 19 * m + 1) for m in range(1, k // 19 + 1))
        for _, _, _, w in words:
            u = GR.from_increments(w)
            good &= len(u) == k and max(w) == 2 and min(w) == 0 and PF.pattern(u) == tuple(w)
            good &= GR.zdim(k, GR.tight_graph(u)[1]) == 3
            if k <= 25:
                good &= z_dim(u) == 3
        ok &= good; rows.append({"k": k, "rays": len(ws)})
    res["counts"] = rows
    checks["lower_bound"] = report(f"padded corners: distinct height-two rays, count sum 2^m (k - 19m + 1), 19 <= k <= {a.kmax}",
                                   ok, f"{rows[-1] if rows else None} ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "height_two_rays")


if __name__ == "__main__":
    main()
