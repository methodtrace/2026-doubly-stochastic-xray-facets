#!/usr/bin/env python3
"""thm:F (facet criterion), thm:empty (empty anti-diagonals), thm:symdual (symmetric dual potential), prop:raybip
and cor:facetbip, against the certified facet test.

  1. Random potentials with a tight permutation, 4 <= k <= --kmax (default 14): integer potentials with small values
     that happen to have a tight permutation, and integer optimal duals (exact Hungarian method) of random
     assignment problems with Hankel weights, kept when they have a tight permutation. For each: the criterion of thm:F (every
     anti-diagonal carries an allowed edge, and dim Z*(u) - bip(u) = 3) agrees with the certified facet test
     (`facet_test.classify`); under (i) the exact face dimension equals (2k - 1) - (dim Z*(u) - bip(u)); dim Z*(u) -
     bip(u) >= 2; dim Y*(u) = bip(u) (computed separately by exact rank).
  2. Every ray of the fan with 5 <= k <= --rays-kmax (default 9; data/rays): the same comparisons, and prop:raybip
     (a ray with a tight permutation has no bipartite component of its tight graph).
  3. thm:empty(b) on all potentials of 1 and 2: whenever some anti-diagonal carries no allowed edge, the exact face
     dimension is below 2k - 4; so every facet among these potentials has an allowed edge on every anti-diagonal.
     cor:facetbip: a facet without total support has bip(u) >= 1 and dim Z*(u) >= 4.
  4. thm:symdual on examples: the potentials of ex:noTS and rem:samefacet (k = 19), two further potentials at k = 17
     (data/counterexamples), and every potential of 1 with an allowed edge on every anti-diagonal. The potential w is
     built as in the proof (strong components of the switching digraph, exact rationals), and every conclusion is
     verified directly: same sum, w(i) + w(j) >= v(i + j) with equality exactly on the allowed graph, w box w = v,
     E(w) = E*(u), total support of w, u - w in Y*(u).
All ranks are exact (integer elimination); the facet test carries its own certificates.
Default: --kmax 14, --samples 40 per k and generator, rays k <= 9 (about 15 s); --full: 200 samples per k and
generator, k <= 16 (a few minutes).
"""
import argparse
import os
import random

from common import DATA, Timer, finish, load_rays, out_dir, read_json, report

from dsxray.facet_criterion import allowed_graph, check_symmetric_dual, components, criterion, y_dim
from dsxray.facet_test import classify, perfect_matching, tight_edges


def random_values(k, rng):
    H = rng.choice([1, 1, 2, 2, 3, 4, 6])
    while True:
        u = [rng.randint(0, H) for _ in range(k)]
        if len(set(u)) > 1 and perfect_matching(k, tight_edges(u)) is not None:
            return u


def assignment_duals(w):
    """Integer optimal duals (alpha, beta) of the assignment problem max sum_i w[i][sigma(i)]: alpha_i + beta_j >=
    w[i][j] for all i, j, with equality on an optimal permutation (Hungarian method with potentials, exact integer
    arithmetic)."""
    n = len(w); INF = None
    cost = [[-w[i][j] for j in range(n)] for i in range(n)]
    u = [0] * (n + 1); v = [0] * (n + 1); p = [0] * (n + 1); way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i; j0 = 0; minv = [INF] * (n + 1); used = [False] * (n + 1)
        while True:
            used[j0] = True; i0 = p[j0]; delta = INF; j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if minv[j] is None or cur < minv[j]:
                        minv[j] = cur; way[j] = j0
                    if delta is None or minv[j] < delta:
                        delta = minv[j]; j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta; v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]; p[j0] = p[j1]; j0 = j1
    alpha = [-u[i] for i in range(1, n + 1)]; beta = [-v[j] for j in range(1, n + 1)]
    assert all(alpha[i] + beta[j] >= w[i][j] for i in range(n) for j in range(n))
    return alpha, beta


def random_dual(k, rng):
    """alpha + beta for integer optimal duals (alpha, beta) of max sum_i v(i + sigma(i)) with random integer v
    (Hankel weights); the potential u = alpha + beta is kept when it has a tight permutation."""
    while True:
        R = rng.choice([3, 6, 10])
        mode = rng.random()
        if mode < 0.4:
            v = [rng.randint(0, R) for _ in range(2 * k - 1)]
        elif mode < 0.7:
            v = [-((R * (t - k) ** 2) // (20 * k)) + rng.randint(0, R) for t in range(2 * k - 1)]
        else:
            v = [(R * (t - k) ** 2) // (20 * k) + rng.randint(0, R) for t in range(2 * k - 1)]
        alpha, beta = assignment_duals([[v[i + j] for j in range(k)] for i in range(k)])
        u = [alpha[i] + beta[i] for i in range(k)]
        if len(set(u)) > 1 and perfect_matching(k, tight_edges(u)) is not None:
            return u


def compare(u, stats, bad, exact_dim=True):
    k = len(u); c = criterion(u, face_dim=exact_dim); f = classify(list(u))["facet"]
    E, A, _ = allowed_graph(u)
    key = "facet" if f else "nonfacet"
    stats[key] = stats.get(key, 0) + 1
    if c["criterion_facet"] != f:
        bad.append(("criterion", u))
    if A is None:
        return c, f
    if c["zdim_star"] - c["bip"] < 2:
        bad.append(("lower_bound_2", u))
    if y_dim(k, A) != c["bip"]:
        bad.append(("dimY_equals_bip", u))
    if exact_dim:
        if c["cover"] and c["face_dim"] != c["formula_dim"]:
            bad.append(("dimension_formula", u))
        if not c["cover"] and c["face_dim"] >= 2 * k - 4:       # thm:empty(b)
            bad.append(("empty_antidiagonal", u))
    elif not c["cover"] and f:
        bad.append(("empty_antidiagonal", u))
    if f and not c["cover"]:
        bad.append(("facet_without_cover", u))
    if f and not c["total_support"]:
        stats["facet_without_total_support"] = stats.get("facet_without_total_support", 0) + 1
        if not (c["bip"] >= 1 and c["zdim_star"] >= 4):                 # cor:facetbip
            bad.append(("corollary_5_2", u))
    if not c["total_support"]:
        stats["no_total_support"] = stats.get("no_total_support", 0) + 1
        if c["cover"]:
            stats["no_total_support_with_cover"] = stats.get("no_total_support_with_cover", 0) + 1
    if c["n_empty"] == 1:
        stats["exactly_one_empty_antidiagonal"] = stats.get("exactly_one_empty_antidiagonal", 0) + 1
    return c, f


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/facet_criterion")
    ap.add_argument("--kmax", type=int, default=14)
    ap.add_argument("--samples", type=int, default=40)
    ap.add_argument("--rays-kmax", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 16); a.samples = max(a.samples, 200)
    out = out_dir(a.out); res = {}; checks = {}; rng = random.Random(a.seed)

    # 1. random potentials with a tight permutation
    t = Timer(); bad = []; stats = {}; pots = []
    for k in range(4, a.kmax + 1):
        for gen in (random_values, random_dual):
            for _ in range(a.samples):
                u = gen(k, rng); pots.append(u)
                compare(u, stats, bad)
    res["random"] = {"potentials": len(pots), "stats": stats, "failures": bad[:20]}
    checks["criterion_random"] = report(f"thm:F = certified facet test, dimension formula, dim Y* = bip, thm:empty(b), cor:facetbip: {len(pots)} random potentials with a tight permutation, 4 <= k <= {a.kmax}",
                                        not bad, f"{stats} ({t()} s)")

    # 2. all rays k <= rays-kmax
    t = Timer(); bad = []; stats = {}; prop = 0; nr = 0
    for k in range(5, a.rays_kmax + 1):
        for u in load_rays(k):
            u = list(u); nr += 1
            compare(u, stats, bad)
            E, A, _ = allowed_graph(u)
            if A is not None and any(b for _, b in components(k, E)):
                prop += 1
    res["rays"] = {"rays": nr, "stats": stats, "failures": bad[:20], "proposition_5_1_failures": prop}
    checks["criterion_rays"] = report(f"thm:F = certified facet test on all {nr} rays with 5 <= k <= {a.rays_kmax}; thm:empty(b); prop:raybip",
                                      not bad and prop == 0, f"{stats} ({t()} s)")

    # 3. thm:symdual on examples
    t = Timer(); ce = read_json(os.path.join(DATA, "counterexamples", "potentials.json"))["potentials"]
    named = [(p["name"], p["u"]) for p in ce if p["name"] in ("example_5_1_chi", "remark_5_2_chi_prime", "k17_total_support_bipartite_component", "k17_cover_without_total_support")]
    results = {}; ok = True; napp = 0
    for name, u in named:
        r = check_symmetric_dual(u); results[name] = r
        ok &= r.get("applicable", False) and all(r.values())
    for u in pots:
        r = check_symmetric_dual(u)
        if r.get("applicable"):
            napp += 1; ok &= all(r.values())
    res["symmetric_dual"] = {"examples": results, "random_applicable": napp}
    checks["symmetric_dual"] = report("thm:symdual: w built and every conclusion verified", ok,
                                      f"{len(named)} named examples, {napp} random potentials with an allowed edge on every anti-diagonal ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "facet_criterion")


if __name__ == "__main__":
    main()
