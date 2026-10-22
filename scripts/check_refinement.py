#!/usr/bin/env python3
"""thm:refine (grid refinement), cor:finer (finer grids) and lem:sprtight, lem:blockgraph, lem:flatlevel, on every
facet ray.

For every facet ray u with 3 <= k <= --kmax (normalized pattern; data/rays and data/facets) and m = 2, 3:
  1. Hypotheses: u is a ray whose tight graph is connected and has total support (checked, not assumed).
  2. lem:flatlevel (a flat level): some (a, b) in E(u), a, b >= 1, has v(a + b - 1) = v(a + b); the tight permutation
     fixing 0 used in the proof exists.
  3. lem:blockgraph: the graph Theta(u) on {1..K} is connected.
  4. lem:sprtight: v' = v(ceil(t/m)) and the tight graph of spr_m u is the predicted one.
  5. thm:refine: spr_m u is normalized with the same H, a ray (exact), with connected tight graph and total support,
     defines a facet of hat-P(mK + 1) (certified test), satisfies P_L, and P_R whenever u does; the permutation pi of
     the proof is tight.
  6. cor:finer: ref_m u = padL^(m-1) spr_m u, a ray with connected tight graph and total support, defining a facet
     of hat-P(mk).
  7. The remark after cor:finer: if the first increment of u is nonzero, spr_m u begins with a nonzero increment
     and ends with an increment below H (so it is no end padding).
Default: k <= 8, m = 2, 3 (300 facet rays, about 1 minute); --full: k <= 9 (about 1.5 minutes).
"""
import argparse

from common import Timer, finish, load_facet_rows, load_rays, out_dir, report

from dsxray import padding as PD
from dsxray import refinement as RF
from dsxray.facet_criterion import graph_z_dim
from dsxray.facet_test import classify, perfect_matching, tight_edges


def is_ray(u):
    return graph_z_dim(len(u), tight_edges(u)) == 3


def facet_patterns(k):
    if k <= 4:
        return sorted(PD.small_facet_patterns(k))
    if k <= 8:
        return sorted({PD.pattern(u) for u in load_rays(k) if classify(list(u))["facet"]})
    return sorted({PD.pattern(r["u"]) for r in load_facet_rows(k) if r["facet"]})


def tight_permutation_fixing_zero(u):
    k = len(u); E = tight_edges(u)
    if (0, 0) not in E:
        return None
    E0 = {(i - 1, j - 1) for (i, j) in E if i >= 1 and j >= 1}
    m = perfect_matching(k - 1, E0)
    if m is None:
        return None
    sigma = [0] * k
    for j, i in m.items():
        sigma[i + 1] = j + 1
    return sigma


def check_one(p, m, fails):
    u = PD.from_pattern(p); H = max(p)
    ts, conn = PD.connected_with_total_support(u)
    if not (is_ray(u) and ts and conn):
        fails["hypotheses"] += 1; return
    sigma = tight_permutation_fixing_zero(u)
    if sigma is None or RF.flat_level(u) is None:
        fails["flat_level"] += 1
    if not RF.theta_connected(u):
        fails["theta"] += 1
    if not RF.tight_graph_lemma_holds(u, m):
        fails["tight_graph_lemma"] += 1
    up = RF.spr(u, m); pp = PD.pattern(up)
    ok = tuple(up[i + 1] - up[i] for i in range(len(up) - 1)) == pp and max(pp) == H
    ts2, conn2 = PD.connected_with_total_support(up)
    ok &= is_ray(up) and ts2 and conn2 and classify(up)["facet"] and PD.P_L(pp) and (PD.P_R(pp) or not PD.P_R(p))
    if sigma is not None:
        pi = RF.refined_tight_permutation(sigma, m); Ep = tight_edges(up)
        ok &= sorted(pi) == list(range(len(up))) and all((i, pi[i]) in Ep for i in range(len(up)))
    if not ok:
        fails["theorem"] += 1
    uq = RF.ref(u, m); pq = PD.pattern(uq); q = pp
    for _ in range(m - 1):
        q = PD.pad_left(q)
    ts3, conn3 = PD.connected_with_total_support(uq)
    if not (pq == q and is_ray(uq) and ts3 and conn3 and classify(uq)["facet"]):
        fails["corollary"] += 1
    if p[0] != 0 and not (pp[0] != 0 and pp[-1] < H):
        fails["core_remark"] += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/refinement")
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--m", default="2,3")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 9)
    out = out_dir(a.out); res = {}; checks = {}
    ms = [int(x) for x in a.m.split(",")]
    t = Timer(); per_k = {}
    fails = dict(hypotheses=0, flat_level=0, theta=0, tight_graph_lemma=0, theorem=0, corollary=0, core_remark=0)
    n = 0
    for k in range(3, a.kmax + 1):
        F = facet_patterns(k); per_k[k] = len(F)
        for p in F:
            for m in ms:
                n += 1; check_one(p, m, fails)
        print(f"   k={k}: {len(F)} facet rays ({t()} s)", flush=True)
    res = {"facet_rays_per_k": per_k, "m": ms, "pairs": n, "failures": fails}
    checks["hypotheses"] = report(f"facet rays 3 <= k <= {a.kmax}: ray, connected tight graph, total support", fails["hypotheses"] == 0)
    checks["lemma_5_12"] = report("lem:flatlevel: a flat level (and a tight permutation fixing 0)", fails["flat_level"] == 0)
    checks["lemma_5_11"] = report("lem:blockgraph: Theta(u) connected", fails["theta"] == 0)
    checks["lemma_5_10"] = report(f"lem:sprtight: v' and the tight graph of spr_m u, m = {a.m}", fails["tight_graph_lemma"] == 0)
    checks["theorem_5_6"] = report(f"thm:refine: spr_m u normalized ray, connected, total support, facet, P_L, P_R; pi tight ({n} pairs)", fails["theorem"] == 0)
    checks["corollary_5_5"] = report("cor:finer: ref_m u = padL^(m-1) spr_m u is a facet ray with total support", fails["corollary"] == 0)
    checks["core_remark"] = report("remark after cor:finer: nonzero first increment gives no end padding", fails["core_remark"] == 0, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "refinement")


if __name__ == "__main__":
    main()
