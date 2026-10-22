#!/usr/bin/env python3
"""Statement E-flat: computer-assisted proof for r <= 6 (thm:Eflat6, app:cells) and exhaustive checks on integer
tails.

  1. Cell enumeration (app:cells, tab:cells): the N = 0, 1, 3, 7, 13, 22 forms E_a + E_b - E_c - E_d (r = 1..6), all
     nonempty cells of their sign arrangement with the unit steps fixed, decided by exact Fourier-Motzkin elimination;
     E-flat is evaluated on each cell. Expected nonempty cells 2, 6, 26, 198, 2,298, 45,014 for r = 1..6 with the
     breakdown of tab:cells (Hall fails, K* empty, E-flat holds) and no failure. The cell lists are compared with
     data/cells/cells_rN.json.
  2. Check (i) of app:cells: every integer tail with r parts and sum at most the safeguard bound lies in an enumerated
     cell, and the verdict of the cell equals the verdict computed directly from the integers; with the bounds of the
     paper (30, 30, 30, 30, 28, 21) every enumerated cell for r <= 4 contains such a tail.
  3. Check (ii) of app:cells (FLOAT): the elimination routine agrees with a floating-point linear program on 3,000
     random small systems and on every branch discarded by the enumeration for r <= --lp-rmax.
  4. E-flat on every tail with g_1 + ... + g_r = n for n <= --nmax (33,554,431 tails for n <= 25, of which 8,177,376
     satisfy Hall's condition with K* nonempty; sec:census, tab:evidence), compared with data/summary/eflat.json.
  5. Consistency with the tight graph: for every staircase with a tight permutation (k <= --kmax), the largest
     balanced set computed from the attainment graph equals the largest closed set of starts of the switching
     digraph, and E-flat agrees with total support.
Defaults: r <= 5 (safeguard sums 16, 16, 16, 16, 14), LP cross-check r <= 4, n <= 16, k <= 12.
--full: r <= 6 with sums 30, 30, 30, 30, 28, 21 (about 20 minutes), LP cross-check r <= 5, n <= 25 (about 30
minutes), k <= 15.
"""
import argparse
import os

import gzip
import json

from common import DATA, Timer, finish, out_dir, read_json, report

from dsxray import cell_enumeration as CE
from dsxray.staircases import staircases

FORMS = {1: 0, 2: 1, 3: 3, 4: 7, 5: 13, 6: 22}
EXPECTED = {1: (2, 0, 2, 0), 2: (6, 0, 6, 0), 3: (26, 0, 20, 6), 4: (198, 12, 90, 96), 5: (2298, 586, 434, 1278),
            6: (45014, 24088, 2654, 18272)}
SAFEGUARD_FULL = {1: 30, 2: 30, 3: 30, 4: 30, 5: 28, 6: 21}
SAFEGUARD_QUICK = {1: 16, 2: 16, 3: 16, 4: 16, 5: 14, 6: 12}


def write_cells(path, obj, gz=False):
    """One cell per line: {"unit_steps": [...], "signs": [...], "verdict": ...}; gzip for large r."""
    text = ('{"r": %d,\n "forms": %s,\n "cells": [\n' % (obj["r"], json.dumps(obj["forms"]))
            + ",\n".join("  " + json.dumps(c) for c in obj["cells"]) + "\n ]}\n")
    if gz:
        with gzip.GzipFile(path + ".gz", "wb", mtime=0) as f:
            f.write(text.encode())
    else:
        with open(path, "w") as f:
            f.write(text)


def deposited_cells(r):
    for name in (f"cells_r{r}.json", f"cells_r{r}.json.gz"):
        p = os.path.join(DATA, "cells", name)
        if os.path.exists(p):
            return read_json(p)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/cells")
    ap.add_argument("--rmax", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=16)
    ap.add_argument("--kmax", type=int, default=12)
    ap.add_argument("--lp-rmax", type=int, default=4)
    ap.add_argument("--write-cells", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 6); a.lp_rmax = max(a.lp_rmax, 5); a.nmax = max(a.nmax, 25); a.kmax = max(a.kmax, 15)
    sg = SAFEGUARD_FULL if a.full else SAFEGUARD_QUICK
    out = out_dir(a.out); res = {"cells": {}, "eflat_tails": {}}; checks = {}
    dep_tails = {r["n"]: r["counts"] for r in read_json(os.path.join(DATA, "summary", "eflat.json"))["eflat_tails"]}
    for r in range(1, a.rmax + 1):
        t = Timer(); C = CE.CellEnumeration(r); lp = {"discarded": 0, "lp_feasible": 0}

        def feas(cons, nv, lp=lp, check=r <= a.lp_rmax):
            ok = CE.feasible(cons, nv)
            if check and not ok:
                lp["discarded"] += 1; lp["lp_feasible"] += CE.lp_feasible_float(cons, nv)
            return ok
        C.run(feasibility=feas); s = C.summary()
        nt, realized, miss = C.safeguard(sg[r])
        s.update(safeguard_sum=sg[r], safeguard_tails=nt, safeguard_cells_realized=realized, safeguard_missing=miss)
        if r <= a.lp_rmax:
            s["discarded_branches_lp_check"] = dict(lp)
        res["cells"][r] = s
        exp = (s["cells"], s["hall_fails"], s["kstar_empty"], s["eflat_holds"])
        ok = exp == EXPECTED[r] and s["eflat_fails"] == 0 and miss == 0 and s["forms"] == FORMS[r]
        if r <= a.lp_rmax:
            ok = ok and lp["lp_feasible"] == 0
        if a.full and r <= 4:
            ok = ok and realized == s["cells"]
        dep = deposited_cells(r)
        if dep is not None:
            ok = ok and dep == C.export()
        checks[f"cells_r{r}"] = report(f"r={r} cell enumeration", ok,
                                       f"{s['cells']} cells, E-flat holds on all with K* nonempty; {nt} tails up to sum {sg[r]} in enumerated cells ({t()} s)")
        if a.write_cells:
            write_cells(os.path.join(out_dir(os.path.join(out, "cells")), f"cells_r{r}.json"), C.export(), gz=r >= 6)
    t = Timer(); agree, dis = CE.random_systems_crosscheck(3000)
    res["random_systems_lp_check"] = {"agree": agree, "disagree": dis}
    checks["elimination_vs_lp"] = report("FLOAT cross-check: elimination agrees with an LP on 3,000 random systems and "
                                         f"rejects only LP-infeasible branches for r <= {a.lp_rmax}", dis == 0 and agree == 3000,
                                         f"({t()} s)")
    t = Timer(); fails = 0
    for n in range(1, a.nmax + 1):
        c = CE.eflat_counts(n); res["eflat_tails"][n] = c["counts"]; fails += c["counts"].get("FAIL", 0)
        if c["counts"] != dep_tails[n]:
            fails += 1
    ntails = sum(sum(v.values()) for v in res["eflat_tails"].values())
    nok = sum(v.get("ok", 0) for v in res["eflat_tails"].values())
    res["eflat_tails_total"] = {"tails": ntails, "hall_and_kstar_nonempty": nok}
    if a.nmax == 25:
        fails += (ntails, nok) != (33554431, 8177376)
    checks["eflat_tails"] = report(f"E-flat on all tails with sum <= {a.nmax}", fails == 0,
                                   f"{ntails} tails, {nok} with Hall's condition and K* nonempty ({t()} s)")
    t = Timer(); cnt = {"ts": 0, "nts": 0}
    for k in range(4, a.kmax + 1):
        for g in staircases(k):
            v = CE.largest_closed_start_set_matches(g)
            if v:
                cnt[v] += 1
    res["consistency"] = dict(cnt, k_range=[4, a.kmax])
    checks["consistency"] = report(f"largest balanced set = largest closed start set, k <= {a.kmax}", True, f"{cnt} ({t()} s)")
    res["checks"] = checks
    finish(res, out, "cells")


if __name__ == "__main__":
    main()
