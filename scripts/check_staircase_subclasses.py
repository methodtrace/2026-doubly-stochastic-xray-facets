#!/usr/bin/env python3
"""Infinite classes of staircases where facets have total support (E-flat), checked on finite ranges.

Claims ((a)-(c) are prop:Eflatclasses and (d) is prop:stairbip, proved in the paper; (e) and (f) are further examples
not cited in the paper; the computations check the constructions on finite ranges):
  (a) unimodal tails: every nonempty balanced set K has 0 in A(K) (no Hall assumption needed); a row q with 0 not
      in A(q) has q at least every maximizer of the prefix average and {q, q + 1} in A(q); hence E-flat holds with the
      rigid level 0, and a staircase with a unimodal tail is a facet only if it has total support;
  (b) unique maximizing pairs: a rigid attaining neighbour of the least row of a balanced set certifies E-flat; so
      E-flat holds when every level has a unique maximizing pair, in particular without additive ties, e.g. for
      every ordering of the superincreasing parts 2, 4, ..., 2^r;
  (c) superadditive prefixes (H_p + H_q <= H_{p+q} for p + q <= r - 1): every row sees column 0, the level 0 is
      attained only by (0, 0), and without total support the anti-diagonal 2E_0 carries no allowed edge;
  (d) every staircase with a tight permutation but without total support has a bipartite allowed component:
      with X = {E_q + 1 : q in K*} and N = N(X), y = 1_X - 1_N is a nonzero vector with y(i) + y(j) = 0 on every
      allowed edge;
  (e) the smallest block of A* need not have a rigid partner: tail (4,5,4,1,5,5,5,5,1,1);
  (f) an example outside the exhaustive ranges: g = (2, 2,3,2,3,2,3,2,3,2,3, 4), k = 31, r = 11, tail sum 29.
Checks (exact; balanced sets by enumerating all row subsets, allowed edges by forcing each edge into a perfect
matching, ranks over Q):
  1. (a) all unimodal tails of length <= --unimodal-rmax with parts in {1..4}; verdicts compared with
     cell_enumeration.eflat_tail (matching-based);
  2. (b) all tails with parts in {1, 2, 3} and r <= 6 (subset and matching computations of K* agree, rigid
     certificate => eflat_tail gives 'ok'), all permutations of (2, 4, ..., 2^r) for r <= 7 (no additive ties,
     never a balanced first row without rigid neighbour), and 32 random permutations each for r = 8, 9, 12, 16,
     24, 32, 48, 64;
  3. (c), (d) every staircase with 3 <= k <= --kmax (4,072 staircases and 1,010 with a tight permutation but without
     total support for k <= 12), total support compared with staircases.total_support;
  4. the two k = 17 potentials of data/counterexamples (not cited in the paper) and the potential of ex:noTS
     (k = 19): dim Z(A), bip(A), face dimension, lost edges, compared with facet_criterion.criterion; (e) and (f).
Default: unimodal r <= 7, staircases k <= 10 (about 3 s). --full: r <= 9 (7,896 tails), k <= 12 (about 20 s).
"""
import argparse
import random
from itertools import permutations, product

from common import Timer, finish, out_dir, report

from dsxray import cell_enumeration as CE
from dsxray import facet_criterion as FC
from dsxray import staircase_subclasses as SS
from dsxray import staircases as ST

TO_CE = {"Hall_fails": "nohall", "strict_Hall": "ts", "endpoint_certificate": "ok"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/staircase_subclasses")
    ap.add_argument("--unimodal-rmax", type=int, default=7)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.unimodal_rmax = max(a.unimodal_rmax, 9); a.kmax = max(a.kmax, 12)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); cnt = {}; n = 0; ok = True
    for r in range(1, a.unimodal_rmax + 1):
        for tail in product(range(1, 5), repeat=r):
            if not SS.unimodal(tail):
                continue
            try:
                v = SS.check_unimodal(tail)
            except AssertionError:
                ok = False; continue
            ok &= CE.eflat_tail(tail) == TO_CE[v]
            cnt[v] = cnt.get(v, 0) + 1; n += 1
    res["unimodal"] = {"rmax": a.unimodal_rmax, "tails": n, "verdicts": cnt}
    if a.unimodal_rmax == 9:
        ok &= n == 7896
    checks["unimodal"] = report(f"(a) unimodal tails, r <= {a.unimodal_rmax}, parts <= 4: column 0 in every balanced set",
                                ok, f"{n} tails {cnt} ({t()} s)")

    t = Timer(); cnt = {}; ok = True; n = 0
    for r in range(1, 7):
        for tail in product(range(1, 4), repeat=r):
            v = SS.rigid_first_row(tail, check_subsets=True); n += 1
            cnt[v] = cnt.get(v, 0) + 1
            if v == "rigid_certificate":
                ok &= CE.eflat_tail(tail) == "ok"
    gen = {}; ng = 0
    for r in range(1, 8):
        for tail in permutations([2 ** j for j in range(1, r + 1)]):
            ok &= SS.no_additive_ties(tail)
            v = SS.rigid_first_row(tail, check_subsets=r <= 6); ng += 1
            gen[v] = gen.get(v, 0) + 1
            ok &= v != "first_row_no_rigid_neighbor"
    rng = random.Random(17)
    for r in (8, 9, 12, 16, 24, 32, 48, 64):
        for _ in range(32):
            tail = [2 ** j for j in range(1, r + 1)]; rng.shuffle(tail)
            ok &= SS.no_additive_ties(tail)
            v = SS.rigid_first_row(tail, check_subsets=r <= 9); ng += 1
            gen[v] = gen.get(v, 0) + 1
            ok &= v != "first_row_no_rigid_neighbor"
    res["rigid_neighbour"] = {"small_tails": n, "small_verdicts": cnt, "superincreasing_tails": ng, "superincreasing_verdicts": gen}
    checks["rigid_neighbour"] = report("(b) rigid neighbour of the first balanced row; superincreasing parts have no ties", ok,
                                       f"{n} small tails, {ng} superincreasing tails {gen} ({t()} s)")

    t = Timer(); c = {"staircases": 0, "no_tight_permutation": 0, "total_support": 0, "bipartite_component": 0,
                      "superadditive": 0}; ok = True
    for k in range(3, a.kmax + 1):
        for g in SS.compositions(k):
            c["staircases"] += 1
            try:
                v = SS.bipartite_component_check(g)
            except AssertionError:
                ok = False; continue
            c["superadditive"] += SS.superadditive_prefix(g)
            c["no_tight_permutation" if v is None else v] += 1
            ok &= (v == "total_support") == ST.total_support(tuple(g))
    res["staircases"] = dict(c, k_range=[3, a.kmax])
    if a.kmax == 12:
        ok &= c["staircases"] == 4072 and c["bipartite_component"] == 1010
    checks["bipartite_and_superadditive"] = report(f"(c), (d) all staircases 3 <= k <= {a.kmax}", ok, f"{c} ({t()} s)")

    t = Timer(); ok = True; ex = {}
    for name, (u, want) in SS.EXAMPLES.items():
        got = SS.example_diagnostics(u)
        cr = FC.criterion(u, face_dim=True)
        ok &= got == want and (cr["zdim_star"], cr["bip"], cr["face_dim"], cr["lost_edges"]) == want
        ex[name] = got
    ok &= SS.selector_counterexample() and SS.outside_range_example()
    res["examples"] = ex
    checks["examples"] = report("examples: k = 17, 17, 19 diagnostics; (e) selector counterexample; (f) k = 31 example", ok,
                                f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "staircase_subclasses")


if __name__ == "__main__":
    main()
