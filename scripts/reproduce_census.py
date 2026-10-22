#!/usr/bin/env python3
"""Facets of hat-P(k) by the smallest step family that produces them (k = 5..12), from the classified facet files.

A facet ray is attributed to the smallest r' = min(r, k-1-r) such that it is a signed r-step potential (a step set and
its complement give the same facet with opposite signs); "other" counts the facets produced by no step potential.
For k = 5..8 the classification is recomputed from the deposited rays; for k >= 9 the deposited facet files are read
(reproduce_facet_counts.py rechecks them). Default and --full: k <= 12 (a few seconds).
"""
import argparse

from common import Timer, finish, load_facet_rows, load_rays, out_dir, report

from dsxray.fan_walk import census
from dsxray.facet_test import classify, classify_coordinate

EXPECTED = {
    5: (23, 9, 8, 6, 0, 0, 0, 0), 6: (45, 11, 10, 20, 0, 0, 0, 4), 7: (85, 13, 12, 30, 17, 0, 0, 13),
    8: (187, 15, 14, 42, 58, 0, 0, 58), 9: (426, 17, 16, 56, 88, 54, 0, 195), 10: (1033, 19, 18, 72, 126, 176, 0, 622),
    11: (2633, 21, 20, 90, 174, 276, 150, 1902), 12: (7041, 23, 22, 110, 234, 404, 500, 5748)}
COLUMNS = ("facets", "coordinate", "single", "double", "triple", "four-step", "five-step", "other")


def facet_rays(k):
    if k <= 8:
        return [list(u) for u in load_rays(k) if classify(list(u))["facet"]]
    return [r["u"] for r in load_facet_rows(k) if r["facet"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/census")
    ap.add_argument("--kmax", type=int, default=12)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    kmax = 12 if a.full else a.kmax
    res = {"table": {}}; checks = {}
    for k in range(5, kmax + 1):
        t = Timer()
        F = facet_rays(k); nc = sum(classify_coordinate(k, s)["facet"] for s in range(2 * k - 1))
        c = census(k, F, nc); c["facets"] = len(F) + nc
        row = tuple(c.get(col, 0) for col in COLUMNS)
        res["table"][k] = dict(zip(COLUMNS, row))
        checks[f"k{k}"] = report(f"k={k} census", row == EXPECTED[k], " ".join(map(str, row)) + f" ({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "census")


if __name__ == "__main__":
    main()
