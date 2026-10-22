#!/usr/bin/env python3
"""Distinctness of the non-step families (thm:fam(d)) and the patterns of lem:distinct, exact.

Notation as in thm:fam: K = k-1, L = k-2, c = k-5; u_{a,b}(i) = min(i, a) on [0, n] and a+2 at n+1, n+2
(a >= 2, b >= 1, n = a+b = k-3); omega = e_1 - e_L; nu = e_K - e_0 - e_2. The pattern of a potential u is
p(u) = (d - min d) / (max d - min d) for its increments d (lem:distinct). For 9 <= k <= --kmax:
  1. the patterns of u_{a,b}, omega, -omega, nu, -nu and of the reflections u o rho equal the formulas printed in the
     proof of thm:fam(d), and the reflection of omega is -omega;
  2. the k - 5 potentials u_{a,b}, omega, -omega, nu, -nu and the reflections of u_{a,b}, nu, -nu have 2k - 4 pairwise
     distinct patterns, none of them in {0, 1}^(k-1) (so none is the pattern of a step potential);
  3. independently of lem:distinct, for k <= --normal-kmax: the facet normals v (modulo affine sequences and positive
     scaling, compared through their second differences) of these potentials are pairwise distinct and differ from
     the normals of all 2^k - 4 signed step potentials (both signs of every step set J with 1 <= |J| <= k-2).
Also: family (a) has exactly k - 5 members with pairwise distinct patterns for 6 <= k <= --kmax.
Defaults: --kmax 24, --normal-kmax 14 (seconds). --full: --kmax 40, --normal-kmax 16 (about 10 s).
"""
import argparse
import itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce

from common import Timer, finish, out_dir, report

from dsxray import nonstep_families as NF
from dsxray.facet_test import v_of

H = Fr(1, 2)


def pat(u):
    d = [u[i + 1] - u[i] for i in range(len(u) - 1)]; lo, hi = min(d), max(d)
    return tuple(Fr(x - lo, hi - lo) for x in d)


def refl(u):
    return list(reversed(u))


def printed(k):
    c = k - 5; out = {}
    for a in range(2, k - 3):
        b = k - 3 - a
        out[("u", a)] = (H,) * a + (0,) * b + (1, 0)
        out[("u_refl", a)] = (1, 0) + (1,) * b + (H,) * a
    out["omega"] = (1, 0) + (H,) * c + (0, 1)
    out["-omega"] = (0, 1) + (H,) * c + (1, 0)
    out["nu"] = (1, 0, 1) + (H,) * c + (1,)
    out["-nu"] = (0, 1, 0) + (H,) * c + (0,)
    out["nu_refl"] = (0,) + (H,) * c + (0, 1, 0)
    out["-nu_refl"] = (1,) + (H,) * c + (1, 0, 1)
    return {key: tuple(Fr(x) for x in v) for key, v in out.items()}


def computed(k):
    out = {}
    for a in range(2, k - 3):
        u = NF.family_A_u(a, k - 3 - a)
        out[("u", a)] = u; out[("u_refl", a)] = refl(u)
    w = NF.W_u(k); nu = NF.Vplus_u(k)
    out["omega"] = w; out["-omega"] = [-x for x in w]
    out["nu"] = nu; out["-nu"] = [-x for x in nu]
    out["nu_refl"] = refl(nu); out["-nu_refl"] = refl([-x for x in nu])
    return out


def normal_key(u):
    """Second differences of v = u box u, scaled to a primitive integer vector (positive scaling kept)."""
    v = v_of(list(u)); d2 = [v[t + 1] - 2 * v[t] + v[t - 1] for t in range(1, len(v) - 1)]
    g = reduce(gcd, [abs(x) for x in d2]) or 1
    return tuple(x // g for x in d2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/family_distinctness")
    ap.add_argument("--kmax", type=int, default=24)
    ap.add_argument("--normal-kmax", type=int, default=14)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 40); a.normal_kmax = max(a.normal_kmax, 16)
    res = {}; checks = {}
    t = Timer(); ok_printed = ok_distinct = True; sizes = {}
    for k in range(9, a.kmax + 1):
        C = computed(k); P = printed(k)
        ok_printed &= all(pat(C[key]) == P[key] for key in P)
        ok_printed &= pat(refl(NF.W_u(k))) == pat([-x for x in NF.W_u(k)])
        pats = {pat(u) for u in C.values()}
        sizes[k] = len(pats)
        ok_distinct &= len(pats) == 2 * k - 4 and all(Fr(1, 2) in p for p in pats)
    res["distinct_patterns"] = sizes
    checks["printed_patterns"] = report(f"patterns as printed in the proof of thm:fam(d), 9 <= k <= {a.kmax}",
                                        ok_printed, f"({t()} s)")
    checks["distinct_patterns"] = report("2k - 4 pairwise distinct patterns, each containing 1/2", ok_distinct)
    ok_a = all(len({pat(NF.family_A_u(aa, k - 3 - aa)) for aa in range(2, k - 3)}) == k - 5 for k in range(6, a.kmax + 1))
    checks["family_a_count"] = report(f"family (a): k - 5 members with distinct patterns, 6 <= k <= {a.kmax}", ok_a)

    t = Timer(); ok = True; stats = {}
    for k in range(9, a.normal_kmax + 1):
        fam = {normal_key(u) for u in computed(k).values()}
        steps = set()
        for r in range(1, k - 1):
            for J in itertools.combinations(range(1, k), r):
                u = [sum(i >= j for j in J) for i in range(k)]
                steps.add(normal_key(u)); steps.add(normal_key([-x for x in u]))
        stats[k] = {"family_normals": len(fam), "step_normals": len(steps), "common": len(fam & steps)}
        ok &= len(fam) == 2 * k - 4 and not (fam & steps)
    res["normals"] = stats
    checks["distinct_normals"] = report(f"facet normals pairwise distinct and distinct from all step potentials, "
                                        f"9 <= k <= {a.normal_kmax}", ok, f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "family_distinctness")


if __name__ == "__main__":
    main()
