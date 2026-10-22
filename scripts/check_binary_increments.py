#!/usr/bin/env python3
"""Binary increments: the interval-cancellation certificate for fan rays.

Claim (a ray certificate for step potentials; not stated in the paper): let u have increments d in {0, 1} (not constant). Consecutive tight pairs
(a, t - a), (b, t - b), a < b, on an anti-diagonal t give the balance identity
    sum_{a <= i < b} d(i) = sum_{t-b <= i < t-a} d(i),
a signed sum of increments with zero coefficient sum, valid for every z in Z(u) in place of d. Merging increment
positions whenever such an identity has exactly two classes with nonzero total coefficient (they must then carry
equal values) and repeating until nothing changes: if exactly two classes remain, Z(u) = span{1, iota, u} and u is a
ray. The test is sound at every size; its converse is open, and holds for k <= 12 (exhaustive). The seven-point
word (1,0,1,1,0,0) is a ray certified this way that is not a facet (its allowed graph misses an anti-diagonal).
Checks (exact): every binary non-constant word with 3 <= k <= --kmax points: each balance identity holds for d;
every merge joins positions with equal increments; two final classes => dim Z = 3 (rational elimination); and the
certificate holds for every ray in this range (1,970 rays at k = 12); dim Z compared with facet_test.z_dim for
k <= 9; the seven-point word: ray, certified, not a facet (facet_test.classify), with an empty anti-diagonal of the
allowed graph (facet_criterion.criterion).
Default: k <= 10 (about 1 s). --full: k <= 12 (about 3 s).
"""
import argparse
from itertools import product

from common import Timer, finish, out_dir, report

from dsxray import exact_graphs as GR
from dsxray import facet_criterion as FC
from dsxray import potential_families as PF
from dsxray.facet_test import classify, z_dim


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/binary_increments")
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 12)
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer(); ok = True; rows = []
    for k in range(3, a.kmax + 1):
        words = rays = cert = 0
        for d in product((0, 1), repeat=k - 1):
            if len(set(d)) < 2:
                continue
            try:
                zd, c = PF.binary_cancellation(d)
            except AssertionError:
                ok = False; continue
            words += 1; rays += zd == 3; cert += c
            ok &= (not c or zd == 3) and (zd != 3 or c)
            if k <= 9:
                ok &= z_dim(GR.from_increments(d)) == zd
        rows.append({"k": k, "words": words, "rays": rays, "certified": cert})
    res["sweep"] = rows
    if a.kmax >= 12:
        ok &= rows[12 - 3]["rays"] == 1970
    checks["cancellation"] = report(f"binary words, k <= {a.kmax}: two final classes iff ray", ok, f"{rows[-1]} ({t()} s)")
    t = Timer()
    d = (1, 0, 1, 1, 0, 0); u = GR.from_increments(d)
    zd, c = PF.binary_cancellation(d); cl = classify(u); cr = FC.criterion(u)
    ok = zd == 3 and c and cl["ray"] and not cl["facet"] and cr["tight_permutation"] and not cr["cover"]
    checks["seven_point_word"] = report("(1,0,1,1,0,0): a certified ray that is not a facet", ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "binary_increments")


if __name__ == "__main__":
    main()
