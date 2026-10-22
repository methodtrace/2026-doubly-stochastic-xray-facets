#!/usr/bin/env python3
"""Staircases with a nonincreasing tail (thm:dec and cor:dec): exact checks.

For every composition g = (g_0, ..., g_r) of k with 1 <= r <= k-2 and g_1 >= ... >= g_r (g_0 arbitrary):
  1. the certified facet test: u(g) with sign + is a facet, a ray (dim Z = 3) and its tight graph has total support;
     the involution g -> g* maps it to a facet as well;
  2. the objects of the proof: M(s) = E_floor(s/2) + E_ceil(s/2) (by brute force over block pairs), and the
     permutation that is the identity on block 0 and reverses every block l >= 1 is tight;
  3. distinctness: the patterns p(u_J)(i) = 1{i+1 in J} of the members are pairwise distinct, so the members define
     pairwise distinct facets (lem:distinct);
  4. the count: the members number sum_{n=1}^{k-1} p(n) - 1 (193 at k = 12, 371 at k = 14, 2,085 at k = 20), and
     p(n) >= (n/m^2)^(m-1) for 2 <= m <= 6, n <= 80 (the bound in the proof of cor:dec);
  5. the reversed condition does not suffice: g = (1, 1, 2, 3) (k = 7) is not a facet.
Default: k <= 14 (the range stated in the paper). --full: k <= 20 (about 15 minutes on 2 cores).
"""
import argparse
from collections import Counter

from common import Timer, finish, out_dir, pmap, report

from dsxray import staircases as ST
from dsxray.facet_test import classify, tight_edges
from dsxray.padding import pattern

PAPER_COUNTS = {12: 193, 14: 371, 20: 2085}


def partitions(N):
    p = [1] + [0] * N
    for part in range(1, N + 1):
        for n in range(part, N + 1):
            p[n] += p[n - part]
    return p


def proof_objects_ok(g):
    r = len(g) - 1; E = ST.ends(g); u = ST.u_of(g); k = len(u)
    M = [max(E[p] + E[s - p] for p in range(max(0, s - r), min(s, r) + 1)) for s in range(2 * r + 1)]
    ok = all(M[s] == E[s // 2] + E[(s + 1) // 2] for s in range(2 * r + 1))
    S = ST.starts(g); sigma = []
    for i in range(g[0]):
        sigma.append(i)
    for l in range(1, r + 1):
        sigma.extend(range(E[l], S[l] - 1, -1))
    T = tight_edges(u)
    return ok and sorted(sigma) == list(range(k)) and all((i, sigma[i]) in T for i in range(k))


def _task(g):
    c = classify(ST.u_of(g)); cs = classify(ST.u_of(ST.star(g)))
    return {"facet": c["facet"], "ray": c["ray"], "total_support": bool(c.get("total_support")),
            "star_facet": cs["facet"], "proof": proof_objects_ok(g)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/nonincreasing_tails")
    ap.add_argument("--kmax", type=int, default=14)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 20)
    p = partitions(max(a.kmax, 80))
    res = {"per_k": {}}; checks = {}
    t = Timer(); bad = Counter(); count_ok = True; distinct_ok = True
    for k in range(3, a.kmax + 1):
        gs = list(ST.nonincreasing_tails(k))
        rows = pmap(_task, gs, a.workers, chunksize=16)
        expected = sum(p[n] for n in range(1, k)) - 1
        pats = {pattern(ST.u_of(g)) for g in gs}
        shifted = all(pattern(ST.u_of(g)) == tuple(int(i + 1 in ST.J_of(g)) for i in range(k - 1)) for g in gs)
        c = Counter()
        for rw in rows:
            c["members"] += 1
            c["facet_ray_total_support"] += rw["facet"] and rw["ray"] and rw["total_support"]
            c["star_facet"] += rw["star_facet"]; c["proof_objects"] += rw["proof"]
        res["per_k"][k] = dict(c, expected_count=expected, distinct_patterns=len(pats))
        for key in ("facet_ray_total_support", "star_facet", "proof_objects"):
            bad[key] += c["members"] - c[key]
        count_ok &= c["members"] == expected and (k not in PAPER_COUNTS or expected == PAPER_COUNTS[k])
        distinct_ok &= len(pats) == len(gs) and shifted
    checks["facets"] = report(f"every member is a ray with total support and a facet, k <= {a.kmax}",
                              bad["facet_ray_total_support"] == 0, f"({t()} s)")
    checks["star"] = report("g -> g* maps the members to facets", bad["star_facet"] == 0)
    checks["proof_objects"] = report("M(s) and the tight permutation of the proof", bad["proof_objects"] == 0)
    checks["distinct"] = report("patterns 1{i+1 in J} pairwise distinct (distinct facets)", distinct_ok)
    checks["count"] = report("count sum_{n<k} p(n) - 1 (193 at k = 12, 371 at k = 14, 2,085 at k = 20)", count_ok,
                             str({k: v["members"] for k, v in res["per_k"].items()}))
    bound_ok = all(p[n] >= (n / m ** 2) ** (m - 1) for m in range(2, 7) for n in range(1, 81))
    bound_ok &= all(sum(p[n] for n in range(1, k)) - 1 >= p[k - 1] for k in range(3, 81))
    res["partitions_1_to_19"] = p[1:20]
    checks["partition_bound"] = report("p(n) >= (n/m^2)^(m-1) (m <= 6, n <= 80) and the sum is >= p(k-1)", bound_ok)
    rev = classify(ST.u_of((1, 1, 2, 3)))
    res["reversed_example"] = {"g": [1, 1, 2, 3], "facet": rev["facet"], "three_step_rule": ST.triple_rule((1, 1, 2, 3))}
    checks["reversed_example"] = report("g = (1, 1, 2, 3) (nondecreasing tail, k = 7) is not a facet",
                                        not rev["facet"] and not ST.triple_rule((1, 1, 2, 3)))
    res["checks"] = checks
    finish(res, out_dir(a.out), "nonincreasing_tails")


if __name__ == "__main__":
    main()
