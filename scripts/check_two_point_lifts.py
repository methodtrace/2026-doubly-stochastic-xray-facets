#!/usr/bin/env python3
"""Two-point block lifts: all two-letter extensions of a facet ray that keep its tight pairs.

Rule (sec:s_census-lifts; verified by exact computation, not proved in the paper): let w be a normalized non-affine
ray on n >= 3 points with connected tight graph and total support, w(0) = 0, v_w = w box w, s0 = w(1),
A = min_{s <= n-3} (w(s+2) - v_w(s)), B = min_{s <= n-2} (w(s+1) - v_w(s)), and for real p, q let u = (-p - q, -q, w(0), ..., w(n-1)) (increments
(p, q, increments of w); negative p, q allowed, normalized afterwards). All shifted tight pairs of w stay tight for u
iff p + q <= A and q <= B; then u is a ray iff (p, q) lies on at least two of the lines p = q, p = s0, p + q = A,
q = B, and every such ray has connected tight graph and total support and defines a facet meeting the left block
face in the facet of w. So there are at most six two-letter extensions: the feasible points among (s0, s0),
(A/2, A/2), (B, B), (s0, A - s0), (s0, B), (A - B, B). With the mirror rule they generate 962 of the 2,808 cores at
k = 12 (536 as left lifts, 536 as right lifts, 110 both; 880 of them not step potentials), exactly the cores whose
block lifts with a block of size 2 are recorded in data/cores/block_lifts_k12.json; the other 1,846 cores are not
two-point lifts.
Checks (exact, Fractions; tight graphs, forced-edge matchings and ranks from exact_graphs.py):
  1. for every facet ray w on n = 5..10 points (data/facets for n >= 8; for n < 8 the rays of data/rays with total
     support): w is a connected ray with total support; each feasible point gives a connected ray with total
     support keeping the shifted pairs; for n + 2 <= 10 each child is a facet ray of data/facets (or a ray with total
     support of data/rays); at n = 8: 172 parents with 1, 2, 3 points for 43, 122, 7 of them and 308 distinct
     children; at n = 10: 1,014 parents;
  2. for the first --grid-parents parents of each n: on a half-integer grid around the lines, the shifted pairs stay
     tight iff p + q <= A and q <= B, and u is a ray iff (p, q) is on at least two lines;
  3. at k = 12: the children of the n = 10 facet rays among the cores (left lifts), their mirrors (right lifts), the
     counts 536, 536, 110, 962, 880, the union equal to the cores with a recorded block lift of size 2, and 1,846
     cores left over; with p, q >= 0 only, 360 left, 360 right and 676 cores.
Default: grid on 3 parents per n (about 20 s). --full: 8 parents per n (about 20 s).
"""
import argparse
import os
from collections import Counter

from common import DATA, Timer, finish, load_facet_rows, load_rays, out_dir, read_json, report

from dsxray import exact_graphs as GR
from dsxray import potential_families as PF


def facet_parents(n):
    if n >= 8:
        rows = [r["u"] for r in load_facet_rows(n) if r["facet"]]
    else:
        rows = [list(u) for u in load_rays(n) if PF.ray_ts_connected(list(u))[1]]
    return sorted({PF.pattern(u) for u in rows})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/two_point_lifts")
    ap.add_argument("--grid-parents", type=int, default=3)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.grid_parents = max(a.grid_parents, 8)
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer(); ok = True; rows = []; children_by_n = {}; grid = 0
    for n in range(5, 11):
        parents = facet_parents(n)
        known = facet_parents(n + 2) if n + 2 <= 10 else None
        known = set(known) if known else None
        dist = Counter(); ch = set()
        for i, wd in enumerate(parents):
            w = GR.from_increments(wd)
            ray, ts, conn, _ = PF.ray_ts_connected(w)
            ok &= ray and ts and conn
            try:
                kids, tested = PF.two_point_lifts(w, grid=i < a.grid_parents)
            except AssertionError:
                ok = False; continue
            grid += tested; dist[len(kids)] += 1; ch.update(kids)
            if known is not None:
                ok &= all(c in known for c in kids)
        children_by_n[n] = ch
        rows.append({"parent_points": n, "parents": len(parents), "points_per_parent": dict(sorted(dist.items())),
                     "distinct_children": len(ch)})
        if n == 10:
            ok &= len(parents) == 1014
        if n == 8:
            ok &= len(parents) == 172 and dict(dist) == {1: 43, 2: 122, 3: 7} and len(ch) == 308
    res["lifts"] = rows; res["grid_points_tested"] = grid
    checks["lifts"] = report("two-point lifts of all facet rays on 5..10 points: rays with total support; grid: (L) and the two-line rule",
                             ok, f"{grid} grid points ({t()} s)")

    t = Timer()
    cores = {tuple(c) for c in read_json(os.path.join(DATA, "cores", "cores_k12.json"))["cores"]}
    left = children_by_n[10] & cores
    right = {PF.mirror(p) for p in left}
    union = left | right
    block = {tuple(r["pattern"]) for r in read_json(os.path.join(DATA, "cores", "block_lifts_k12.json"))["cores"]
             if any(b == 2 for _, b in r["block_lifts"])}
    nonstep = sum(max(p) > 1 for p in union)
    counts = {"cores": len(cores), "left": len(left), "right": len(right), "both": len(left & right), "union": len(union),
              "not_step": nonstep, "other_cores": len(cores - union)}
    ok = (counts == {"cores": 2808, "left": 536, "right": 536, "both": 110, "union": 962, "not_step": 880, "other_cores": 1846}
          and right <= cores and union == block)
    res["cores_k12"] = counts
    checks["cores_k12"] = report("k = 12: two-point lifts give 962 of the 2,808 cores (536 left, 536 right, 110 both)", ok,
                                 f"{counts} ({t()} s)")
    t = Timer(); nonneg = set()
    for wd in facet_parents(10):
        w = GR.from_increments(wd)
        for p, q in PF.lift_candidates(w)[3]:
            if p >= 0 and q >= 0:
                nonneg.add(PF.pattern([-p - q, -q] + list(w)))
    ln = nonneg & cores; rn = {PF.mirror(x) for x in ln}
    res["nonnegative_only"] = {"left": len(ln), "right": len(rn), "union": len(ln | rn)}
    checks["nonnegative_only"] = report("with p, q >= 0 only: 360 left, 360 right, 676 cores", (len(ln), len(rn), len(ln | rn)) == (360, 360, 676),
                                        f"{res['nonnegative_only']} ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "two_point_lifts")


if __name__ == "__main__":
    main()
