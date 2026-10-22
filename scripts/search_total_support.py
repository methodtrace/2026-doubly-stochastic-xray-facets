#!/usr/bin/env python3
"""Exhaustive search over integer potentials for failures of total support (C program tools/c/search_total_support.c).

  (a) Values: all potentials u in {0..H}^k with u(0) = 0. Among those with a tight permutation and without total
      support, the number of anti-diagonals that carry no allowed edge (the histogram and its minimum), and every
      potential with an allowed edge on every anti-diagonal ("cover without total support"), with its ray and facet
      flags. Also every potential with total support whose tight graph has a bipartite component.
  (b) Increments: all normalized potentials of height at most H (increments in {0..H}^(k-1), minimum 0, gcd 1):
      the rays with a tight permutation, those without total support, and every facet ray without total support.
Every potential the C program reports is re-verified in Python in exact arithmetic (facet_criterion.criterion and
the certified facet test). The C program computes ranks modulo two primes: a reported ray or facet is certified,
and a ray or facet could be missed only if both modular ranks fell below the rational rank.
Expected outcomes (exact, from the same program): for (k, H) = (10, 4), (8, 7), (12, 3) no potential has an
allowed edge on every anti-diagonal without total support (minimum 3, 3, 2 empty anti-diagonals); for (19, 1) four
potentials do, and all four are facet rays (the potential of ex:noTS, the potential chi' of rem:samefacet, and two
more); in (b) no facet ray of height at most 2 lacks total support for k <= 18 (rem:samefacet).
Default: (a) for (10, 4), (8, 7), (12, 3), (19, 1) and (b) for k <= 12, H = 2 (about 15 s); --full: (a) also for
(13, 2) and (16, 2) (minimum 2 and 1 empty anti-diagonals, no potential with an allowed edge on every anti-diagonal
without total support) and (b) for k <= 16 (a few minutes). The range k <= 18 of rem:samefacet is run by
check_facets_without_total_support.py --full (or --inc-kmax 18 here; about 28 core-minutes, --workers splits each k by
the last increment).
"""
import argparse
import json
import os
import subprocess
from multiprocessing import Pool

from common import Timer, build_c, finish, out_dir, report

from dsxray.facet_criterion import criterion
from dsxray.facet_test import classify

VALUES_DEFAULT = [(10, 4), (8, 7), (12, 3), (19, 1)]
EXPECTED_VALUES = {(10, 4): {"cover_without_total_support": 0, "min_empty_antidiagonals": 3},
                   (8, 7): {"cover_without_total_support": 0, "min_empty_antidiagonals": 3},
                   (12, 3): {"cover_without_total_support": 0, "min_empty_antidiagonals": 2},
                   (19, 1): {"cover_without_total_support": 4, "min_empty_antidiagonals": 0},
                   (13, 2): {"cover_without_total_support": 0, "min_empty_antidiagonals": 2},
                   (16, 2): {"cover_without_total_support": 0, "min_empty_antidiagonals": 1}}


def _run(args):
    p = subprocess.run([str(x) for x in args], check=True, capture_output=True, text=True)
    lines = p.stdout.splitlines()
    return json.loads(lines[-1]), lines[:-1]


def run_values(exe, k, H):
    summary, lines = _run([exe, "V", k, H])
    cover = []; tsbip = []
    for line in lines:
        if line.startswith("COVER"):
            a, b = line[5:].split("|"); f = b.split()
            cover.append({"u": list(map(int, a.split())), "ray": f[1] == "1", "facet": f[3] == "1"})
        elif line.startswith("TSBIP"):
            tsbip.append(list(map(int, line[5:].split())))
    return summary, cover, tsbip


def run_increments(exe, k, H, workers=1):
    """Search (b) for one k, split into H + 1 runs by the last increment."""
    jobs = [[exe, "I", k, H, last] for last in range(H + 1)]
    if workers > 1:
        with Pool(min(workers, len(jobs))) as pool:
            outs = pool.map(_run, jobs)
    else:
        outs = [_run(j) for j in jobs]
    tot = {"k": k, "H": H, "normalized_potentials": 0, "rays_with_tight_permutation": 0, "rays_without_total_support": 0,
           "facet_rays_without_total_support": 0}
    found = []
    for summary, lines in outs:
        for key in list(tot)[2:]:
            tot[key] += summary[key]
        found += [list(map(int, line[5:].split())) for line in lines if line.startswith("FACET")]
    return tot, found


def verify_cover(entry):
    c = criterion(entry["u"])
    f = classify(entry["u"])
    return (c["tight_permutation"] and not c["total_support"] and c["cover"] and (c["zdim"] == 3) == entry["ray"]
            and f["facet"] == entry["facet"] == c["criterion_facet"])


def verify_facet_without_ts(u):
    c = criterion(u); f = classify(u)
    return c["zdim"] == 3 and not c["total_support"] and f["facet"] and c["criterion_facet"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/search_total_support")
    ap.add_argument("--values", default=",".join(f"{k}:{H}" for k, H in VALUES_DEFAULT),
                    help="comma-separated k:H for search (a)")
    ap.add_argument("--inc-kmax", type=int, default=12)
    ap.add_argument("--inc-H", type=int, default=2)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.inc_kmax = max(a.inc_kmax, 16)
        if a.values == ap.get_default("values"):
            a.values += ",13:2,16:2"
    out = out_dir(a.out); res = {"values": {}, "increments": {}}; checks = {}
    exe = build_c("search_total_support", os.path.join(out, "bin"))

    t = Timer(); ok = True
    for item in a.values.split(","):
        k, H = map(int, item.split(":"))
        summary, cover, tsbip = run_values(exe, k, H)
        good = all(verify_cover(e) for e in cover)
        exp = EXPECTED_VALUES.get((k, H))
        if exp:
            good &= all(summary[key] == val for key, val in exp.items())
        ok &= good
        res["values"][f"{k}:{H}"] = {"summary": summary, "cover_without_total_support": cover,
                                     "total_support_with_bipartite_component": tsbip, "verified": good}
        print(f"   values k={k} H={H}: {summary['tight_permutation']} with a tight permutation, {summary['no_total_support']} without total support, "
              f"min empty anti-diagonals {summary['min_empty_antidiagonals']}, cover without total support {len(cover)} "
              f"({sum(e['facet'] for e in cover)} facets), total support with bipartite component {len(tsbip)}", flush=True)
    checks["values"] = report("(a) potentials with bounded values: outcomes as expected, every reported potential re-verified exactly", ok, f"({t()} s)")

    t = Timer(); ok = True; found_all = []
    for k in range(5, a.inc_kmax + 1):
        tot, found = run_increments(exe, k, a.inc_H, a.workers)
        ok &= all(verify_facet_without_ts(u) for u in found)
        found_all += found
        res["increments"][k] = {"summary": tot, "facet_rays_without_total_support": found}
        print(f"   increments k={k} H={a.inc_H}: {tot['rays_with_tight_permutation']} rays with a tight permutation, "
              f"{tot['rays_without_total_support']} without total support, {len(found)} facet rays without total support ({t()} s)", flush=True)
    if a.inc_H <= 2 and a.inc_kmax <= 18:
        ok &= not found_all
    checks["increments"] = report(f"(b) facet rays of height <= {a.inc_H} without total support, 5 <= k <= {a.inc_kmax}: "
                                  + ("none" if not found_all else f"{len(found_all)} found, re-verified"), ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "search_total_support")


if __name__ == "__main__":
    main()
