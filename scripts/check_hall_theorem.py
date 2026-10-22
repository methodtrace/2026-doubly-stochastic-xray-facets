#!/usr/bin/env python3
"""thm:H (Hall sets of block starts): computational checks of the statement and of every step of its proof.

  1. On all staircases with 3 <= k <= --kmax: total support fails iff some nonempty set of starts S_1..S_r has at most
     |X| tight neighbors, decided by the matching criterion and compared with the certified total-support test; the
     tight neighbors of the starts are block ends (262,108 staircases for k <= 18).
  2. For k <= --brute-kmax the matching criterion agrees with brute force over all sets of starts.
  3. For k <= --proof-kmax every inequality used in the proof is checked on the actual tight graph for every
     down-closed row set (3,012,175 pairs for k <= 13).
  4. Minimal Hall sets of the staircases without total support, k <= --sizes-kmax (sizes 1, 2, 3 in 5,272, 1,949
     and 157 of the 7,378 cases with k <= 14).
Defaults: 1 for k <= 13, 2 for k <= 10, 3 for k <= 10, 4 for k <= 12. --full: 18, 12, 13, 14 (about 1.5 hours).
"""
import argparse
from collections import Counter

from common import Timer, finish, out_dir, pmap, report

from dsxray import hall_sets as HS
from dsxray.staircases import staircases

HALL_COUNTS = {3: 2, 4: 6, 5: 14, 6: 30, 7: 62, 8: 126, 9: 254, 10: 510, 11: 1022, 12: 2046, 13: 4094, 14: 8190,
               15: 16382, 16: 32766, 17: 65534, 18: 131070}
NO_TS = {7: 3, 8: 12, 9: 40, 10: 118, 11: 312, 12: 776, 13: 1849, 14: 4268, 15: 9579, 16: 21086, 17: 45607, 18: 97378}


def _thm(g):
    return HS.check_theorem(g)


def _proof(g):
    return HS.proof_step_checks(g)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/hall_theorem")
    ap.add_argument("--kmax", type=int, default=13)
    ap.add_argument("--brute-kmax", type=int, default=10)
    ap.add_argument("--proof-kmax", type=int, default=10)
    ap.add_argument("--sizes-kmax", type=int, default=12)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax, a.brute_kmax, a.proof_kmax, a.sizes_kmax = max(a.kmax, 18), max(a.brute_kmax, 12), max(a.proof_kmax, 13), max(a.sizes_kmax, 14)
    res = {"theorem": {}}; checks = {}; total = 0
    for k in range(3, a.kmax + 1):
        t = Timer(); out = pmap(_thm, list(staircases(k)), a.workers, chunksize=256)
        n = len(out); nts = sum(not ts for ts, w, e in out); exc = sum(ts == w for ts, w, e in out); ne = sum(not e for _, _, e in out)
        total += n
        res["theorem"][k] = {"staircases": n, "no_total_support": nts, "exceptions": exc, "start_neighbors_not_ends": ne}
        checks[f"theorem_k{k}"] = report(f"k={k} thm:H", exc == 0 and ne == 0 and n == HALL_COUNTS[k]
                                         and nts == NO_TS.get(k, 0), f"{n} staircases, {nts} without total support ({t()} s)")
    res["theorem_total"] = total
    if a.kmax == 18:
        checks["total_262108"] = report("262,108 staircases for k <= 18", total == 262108)
    t = Timer(); C = Counter()
    for k in range(3, a.brute_kmax + 1):
        C.update(HS.theorem_counts(k, brute=True))
    res["brute_force"] = dict(C, k_range=[3, a.brute_kmax])
    checks["brute_force"] = report(f"matching criterion equals brute force over all sets of starts, k <= {a.brute_kmax}",
                                   C["brute_force_disagreements"] == 0 and C["exceptions"] == 0, f"({t()} s)")
    t = Timer(); npairs = 0
    for k in range(3, a.proof_kmax + 1):
        npairs += sum(pmap(_proof, list(staircases(k)), a.workers, chunksize=64))
    res["proof_steps"] = {"k_range": [3, a.proof_kmax], "staircase_row_set_pairs": npairs}
    checks["proof_steps"] = report(f"proof inequalities on every down-closed row set, k <= {a.proof_kmax}", True,
                                   f"{npairs} pairs ({t()} s)")
    if a.proof_kmax == 13:
        checks["proof_pairs_3012175"] = report("3,012,175 pairs for k <= 13", npairs == 3012175)
    t = Timer(); n, sizes = HS.minimal_size_counts(3, a.sizes_kmax)
    res["minimal_hall_sets"] = {"k_range": [3, a.sizes_kmax], "without_total_support": n, "sizes": sizes}
    checks["minimal_sizes"] = report(f"minimal Hall sets, k <= {a.sizes_kmax}", None not in sizes and "None" not in sizes,
                                     f"{n} staircases, sizes {sizes} ({t()} s)")
    if a.sizes_kmax == 14:
        checks["sizes_k14"] = report("sizes 5272, 1949, 157 of 7378", n == 7378 and sizes == {"1": 5272, "2": 1949, "3": 157})
    res["checks"] = checks
    finish(res, out_dir(a.out), "hall_theorem")


if __name__ == "__main__":
    main()
