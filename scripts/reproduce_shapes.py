#!/usr/bin/env python3
"""Shapes of the facet rays: hinges, heights, largest height (k = 4..12).

Every ray is, after subtracting a linear function and scaling, a potential whose increments are nonnegative integers
with minimum 0 and gcd 1; its height is the largest increment (u and -u have the same height). Checks:
  - the convex and the concave facet rays are exactly the 2(k-2) hinges (i - c)^+ and their negatives;
  - the table of facet rays by height (hinges, other height 1, height 2, height 3, height >= 4) for k = 6..12;
  - the largest height 1, 1, 2, 3, 4, 5, 7, 9, 12 for k = 4..12.
Default and --full: k <= 12 (a few seconds).
"""
import argparse
from collections import Counter

from common import Timer, finish, load_facet_rows, load_rays, out_dir, report

from dsxray.facet_test import classify
from dsxray.padding import increments, pattern
from dsxray.small_k_fan import line_enumeration, normalize_ray

TABLE = {6: (34, 8, 22, 4, 0, 0), 7: (72, 10, 49, 12, 1, 0), 8: (172, 12, 102, 48, 8, 2), 9: (409, 14, 200, 143, 39, 13),
         10: (1014, 16, 376, 410, 138, 74), 11: (2612, 18, 692, 1089, 490, 323), 12: (7018, 20, 1250, 2780, 1572, 1396)}
MAX_HEIGHT = {4: 1, 5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 10: 7, 11: 9, 12: 12}


def facet_rays(k):
    if k <= 4:
        rays = {normalize_ray(u) for u in line_enumeration(k)[0]}
        return [list(u) for u in rays if classify(list(u))["facet"]]
    if k <= 8:
        return [list(u) for u in load_rays(k) if classify(list(u))["facet"]]
    return [r["u"] for r in load_facet_rows(k) if r["facet"]]


def hinge_kinks(F, k):
    conv, conc = [], []
    for u in F:
        d = increments(u)
        if all(d[i] <= d[i + 1] for i in range(k - 2)):
            conv.append(d)
        elif all(d[i] >= d[i + 1] for i in range(k - 2)):
            conc.append(d)
    ok = True
    for D in (conv, conc):
        ok &= len(D) == k - 2 and all(len(set(d)) == 2 for d in D)
        if ok:
            ok &= sorted(next(i for i in range(k - 2) if d[i] != d[i + 1]) + 1 for d in D) == list(range(1, k - 1))
    return ok, len(conv) + len(conc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/shapes")
    ap.add_argument("--kmax", type=int, default=12)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    kmax = 12 if a.full else a.kmax
    res = {"table": {}, "height_distribution": {}, "largest_height": {}}; checks = {}
    for k in range(4, kmax + 1):
        t = Timer(); F = facet_rays(k)
        ok, nh = hinge_kinks(F, k)
        checks[f"k{k}_hinges"] = report(f"k={k} convex and concave facet rays are the 2(k-2) hinges", ok)
        H = Counter(max(pattern(u)) for u in F)
        res["height_distribution"][k] = dict(sorted(H.items())); res["largest_height"][k] = max(H)
        checks[f"k{k}_largest_height"] = report(f"k={k} largest height", max(H) == MAX_HEIGHT[k], str(max(H)))
        if k >= 6:
            row = (len(F), nh, H[1] - nh, H[2], H[3], sum(v for h, v in H.items() if h >= 4))
            res["table"][k] = dict(zip(("ray_facets", "hinges", "other_height_1", "height_2", "height_3", "height_4_or_more"), row))
            checks[f"k{k}_table"] = report(f"k={k} heights", row == TABLE[k], " ".join(map(str, row)) + f" ({t()} s)")
    if 12 in res["table"]:
        row = res["table"][12]; h2 = row["height_2"] + row["height_3"] + row["height_4_or_more"]
        checks["k12_repeated_thresholds"] = report("k = 12: 5,748 of the 7,018 ray facets have height >= 2 (rem:multiplicity); "
                                                   "height 2 accounts for 2,780", h2 == 5748 and row["height_2"] == 2780)
    res["checks"] = checks
    finish(res, out_dir(a.out), "shapes")


if __name__ == "__main__":
    main()
