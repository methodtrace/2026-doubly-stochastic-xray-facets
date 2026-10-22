#!/usr/bin/env python3
"""thm:T3 (three steps; cor:T3, tab:evidence): u(g), g = (g0, g1, g2, g3), is a facet unless g1, g2, g3 are
pairwise distinct and g1 < g3.

  1. Every step of the proof by brute force on all three-step staircases with 5 <= k <= --proof-kmax: closed-form
     corner-triangle sizes, the neighbor sets of the starts in the failing cases, and in the facet cases the explicit
     tight involution, the strong connectivity of its switching digraph and the connectivity of both step classes.
  2. The statement against the certified facet test on all three-step staircases with 5 <= k <= --kmax
     (91,389 potentials for k <= 40).
Defaults: proof steps k <= 14, statement k <= 22. --full: proof steps k <= 30 (about 1 hour), statement k <= 40
(about 2.5 core-hours).
"""
import argparse
import itertools

from common import Timer, finish, out_dir, pmap, report

from dsxray import staircases as ST
from dsxray.facet_test import classify


def _stmt(g):
    return classify(ST.u_of(g))["facet"] == ST.triple_rule(g)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/three_step")
    ap.add_argument("--proof-kmax", type=int, default=14)
    ap.add_argument("--kmax", type=int, default=22)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.proof_kmax = max(a.proof_kmax, 30); a.kmax = max(a.kmax, 40)
    checks = {}; res = {}
    t = Timer(); pc = ST.three_step_proof_checks(5, a.proof_kmax, theorem_check=False)
    pc = {key: v for key, v in pc.items() if key != "theorem"}
    res["proof_steps"] = {"k_range": [5, a.proof_kmax],
                          **{key: {"checked": v[0], "passed": v[1]} for key, v in pc.items()}}
    checks["proof_steps"] = report(f"proof steps, 5 <= k <= {a.proof_kmax}",
                                   all(v[0] > 0 and v[0] == v[1] for v in pc.values()),
                                   f"{ {key: v[0] for key, v in pc.items()} } ({t()} s)")
    t = Timer(); n = ok = 0
    for k in range(5, a.kmax + 1):
        gs = [(x, y - x, z - y, k - z) for x, y, z in itertools.combinations(range(1, k), 3)]
        out = pmap(_stmt, gs, a.workers, chunksize=128); n += len(out); ok += sum(out)
    res["statement"] = {"k_range": [5, a.kmax], "potentials": n, "agree": ok}
    expect = sum((k - 1) * (k - 2) * (k - 3) // 6 for k in range(5, a.kmax + 1))
    checks["statement"] = report(f"statement on all three-step potentials, k <= {a.kmax}", n == ok == expect,
                                 f"{ok}/{n} ({t()} s)")
    if a.kmax == 40:
        checks["count_91389"] = report("91,389 three-step potentials for k <= 40", n == 91389)
    res["checks"] = checks
    finish(res, out_dir(a.out), "three_step")


if __name__ == "__main__":
    main()
