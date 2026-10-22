#!/usr/bin/env python3
"""The family g_b = G_{1,b} of thm:height (height b on 2b + 2 points) and the quadratic family, as an exact check.

  1. g_b = G_{1,b} (thm:height with a = 1), for 2 <= b <= --bmax (default 10; --full: 40): g_b on {0..2b+1} with
     increments (b, 0^(b-1), b, 1^b) (k = 2b + 2). With w = g_b - iota, the values v(t) and the complete tight graph stated in the proof; the index b
     is tight with every index; the stated permutation is tight and fixes b; total support (computed), connected
     tight graph; ray (exact dim Z = 3); height b; P_L (increment 0 at 1, (0, 1) tight) and P_R (increment b = H at b,
     (b + 1, 2b + 1) tight); I(g_b) defines a facet of hat-P(2b + 2) (certified facet test); the first increment is
     nonzero and the last is below H, so g_b is no end padding (a core).
  2. The first bound of thm:height: for 3 <= k <= 2 --bmax + 3, some facet ray of hat-P(k) has height >= floor((k - 2)/2): single
     steps for k <= 5, and padL^e g_b (k = 2b + 2 + e, e in {0, 1}) otherwise (certified facet test).
  3. The certified facet test on the quadratic family G_ab (its theorem is checked by check_quadratic_family.py):
     the potentials with increments
     (b^a, 0^(b-1), (b-1)a + 1, 1^b) (k = a + 2b + 1; a = 1 is g_b) define facets exactly for 1 <= a <= b + 1, for
     2 <= b <= --gmax (default 4; --full: 8) and 1 <= a <= b + 3; height b^2 at a = b + 1. The outcomes are compared
     with data/families/quadratic_height_family.json (b <= 8).
Default: seconds; --full: b <= 40 (k = 82) and the family with b <= 8, about 2 minutes.
"""
import argparse
import os

from common import DATA, Timer, finish, out_dir, read_json, report, write_json

from dsxray import nonstep_families as NF
from dsxray import padding as PD
from dsxray.facet_criterion import graph_z_dim
from dsxray.facet_test import classify, tight_edges, v_of


def check_g(b):
    p = NF.height_family_pattern(b); u = PD.from_pattern(p); k = len(u); K = k - 1
    w = [u[i] - i for i in range(k)]
    ok = {}
    ok["v"] = v_of(w) == NF.height_family_v(b)
    E = tight_edges(u); ok["tight_graph"] = tight_edges(w) == E == NF.height_family_tight_graph(b)
    ok["universal_index"] = all((b, j) in E for j in range(k))
    s = NF.height_family_permutation(b)
    ok["permutation"] = sorted(s) == list(range(k)) and s[b] == b and all((i, s[i]) in E for i in range(k))
    ts, conn = PD.connected_with_total_support(u); ok["total_support_connected"] = ts and conn
    ok["ray"] = graph_z_dim(k, E) == 3
    ok["height"] = PD.pattern(u) == p and max(p) == b
    ok["P_L"] = p[1] == 0 and (0, 1) in E and PD.P_L(p)
    ok["P_R"] = p[b] == b == max(p) and (b + 1, K) in E and PD.P_R(p)
    ok["facet"] = classify(u)["facet"]
    ok["core"] = p[0] != 0 and p[-1] < max(p)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/unbounded_height")
    ap.add_argument("--bmax", type=int, default=10)
    ap.add_argument("--gmax", type=int, default=4)
    ap.add_argument("--write-data", action="store_true", help="rewrite data/families/quadratic_height_family.json")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.bmax = max(a.bmax, 40); a.gmax = max(a.gmax, 8)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); bad = {}
    for b in range(2, a.bmax + 1):
        r = check_g(b)
        if not all(r.values()):
            bad[b] = [x for x, y in r.items() if not y]
    res["theorem_5_8_failures"] = bad
    checks["theorem_5_8"] = report(f"thm:height (a = 1): g_b, 2 <= b <= {a.bmax} (k <= {2 * a.bmax + 2}): v, tight graph, universal index, tight permutation, total support, ray, height b, P_L, P_R, facet, core",
                                   not bad, f"({t()} s)")

    t = Timer(); bad = []
    for k in range(3, 2 * a.bmax + 4):
        if k <= 5:
            ok = classify([0] + [1] * (k - 1))["facet"] and (k - 2) // 2 <= 1
        else:
            b = (k - 2) // 2; e = k - 2 - 2 * b; q = NF.height_family_pattern(b)
            for _ in range(e):
                q = PD.pad_left(q)
            ok = len(q) == k - 1 and max(q) == b and classify(PD.from_pattern(q))["facet"]
        if not ok:
            bad.append(k)
    checks["height_bound"] = report(f"some facet ray of height >= floor((k-2)/2) for 3 <= k <= {2 * a.bmax + 3}", not bad, f"({t()} s)")

    t = Timer(); rows = []
    for b in range(2, a.gmax + 1):
        for aa in range(1, b + 4):
            p = NF.quadratic_height_pattern(aa, b); u = PD.from_pattern(p)
            rows.append({"a": aa, "b": b, "k": len(u), "height": max(PD.pattern(u)), "facet": classify(u)["facet"]})
    ok = all(r["facet"] == (r["a"] <= r["b"] + 1) for r in rows)
    ok &= all(r["height"] == r["b"] ** 2 for r in rows if r["a"] == r["b"] + 1)
    path = os.path.join(DATA, "families", "quadratic_height_family.json")
    if a.write_data:
        write_json(path, {"description": "Certified facet test on the quadratic family G_ab (the theorem is checked by scripts/check_quadratic_family.py): the normalized potential with increments (b^a, 0^(b-1), (b-1)a+1, 1^b), k = a + 2b + 1, and whether its inequality defines a facet of hat-P(k) (certified facet test). A facet exactly for 1 <= a <= b + 1; a = 1 is the family g_b (Theorem 5.9).",
                          "rows": rows})
    dep = {(r["a"], r["b"]): r for r in read_json(path)["rows"]}
    ok &= all(dep.get((r["a"], r["b"])) == r for r in rows)
    res["quadratic_height_family"] = rows
    checks["quadratic_height_family"] = report(f"quadratic family G_ab, certified facet test: increments (b^a, 0^(b-1), (b-1)a+1, 1^b): facet iff a <= b+1, 2 <= b <= {a.gmax}; equal to the deposited table",
                                               ok, f"{len(rows)} potentials ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "unbounded_height")


if __name__ == "__main__":
    main()
