#!/usr/bin/env python3
"""Statement (S1) for every sequence with r <= 10 boundary steps (a second route to thm:S1cert), by
the exact enumeration of attaining structures (tools/c/attaining_structures.c), and the catalogue of attaining
structures for r <= 6 (data/attaining_structures).

  1. (S1) mode, r = 1..--rmax: a depth-first search over the levels with exact feasibility decisions, each backed by
     a verified integer witness or a verified Gordan multiplier (the C program aborts if a check fails); a leaf
     reached would be a counterexample to (S1). The counts for r = 8, 9 are compared with the deposited values.
     r = 10 is run in chunks: --r10-chunks 6 runs the six chunks obtained by splitting the branches at level 6
     (--r10-only I,J,... selects some of them); about 2.3 core-hours in all (the chunks take between 10 and 45
     minutes); the counts of each chunk are compared with the deposited values.
  2. Catalogue, r = 1..--catmax (default 6): every realizable attaining structure with an integer witness. The counts
     (1, 3, 13, 75, 555, 4,421 structures; 0, 0, 1, 13, 154, 1,602 with several classes; 1, 3, 11, 46, 228, 1,104
     with every q < r having two partners) are checked, the output is compared with data/attaining_structures, and
     the deposited catalogue is rechecked by an independent Python checker (dsxray.attaining_structures): each
     witness realizes exactly its structure, the structures are distinct, (S1) holds on each. Completeness
     cross-check: the structures of all integer sequences with values in {0..V} lie in the catalogue (r <= 4,
     V = 6; r = 5, V = 4).
Default: (S1) for r <= 7 and the catalogue for r <= 6 (about 30 s). --full: (S1) for r <= 9 (r = 9: about 8 minutes on
one core) and --r10-chunks 0 unless given.
"""
import argparse
import json
import os
import subprocess
from multiprocessing import Pool

from common import DATA, Timer, build_c, finish, out_dir, read_jsonl, report, write_jsonl

from dsxray import attaining_structures as AS

EXPECTED_S1 = {8: {"nodes": 134329, "feasibility_tests": 932856, "infeasibility_certificates": 474756, "connectivity_pruned": 32732},
               9: {"nodes": 1129778, "feasibility_tests": 9952183, "infeasibility_certificates": 4843894, "connectivity_pruned": 232236}}
EXPECTED_R10 = {0: (1937011, 20663125, 9862592, 346003), 1: (2570029, 26933214, 12785340, 477725),
                2: (676039, 7677617, 3357700, 124236), 3: (1632575, 18956974, 8970455, 280251),
                4: (1885816, 20838377, 9644943, 337525), 5: (674316, 8162060, 3546068, 116786)}
EXPECTED_CAT = {1: (1, 0, 1), 2: (3, 0, 3), 3: (13, 1, 11), 4: (75, 13, 46), 5: (555, 154, 228), 6: (4421, 1602, 1104)}
BRUTE = {1: 6, 2: 6, 3: 6, 4: 6, 5: 4}


def run(args):
    p = subprocess.run([str(x) for x in args], check=True, capture_output=True, text=True)
    lines = p.stdout.splitlines()
    return json.loads(lines[-1]), lines[:-1]


def cat_path(r):
    return os.path.join(DATA, "attaining_structures", f"structures_r{r}.jsonl" + (".gz" if r >= 6 else ""))


def load_catalogue(r):
    return [(rec["levels"], rec["witness"]) for rec in read_jsonl(cat_path(r))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/attaining_structures")
    ap.add_argument("--rmax", type=int, default=7)
    ap.add_argument("--catmax", type=int, default=6)
    ap.add_argument("--r10-chunks", type=int, default=0)
    ap.add_argument("--r10-only", default="")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--write-data", action="store_true", help="rewrite data/attaining_structures from this run")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 9)
    out = out_dir(a.out); res = {"s1": {}, "catalogue": {}}; checks = {}
    exe = build_c("attaining_structures", os.path.join(out, "bin"))

    ok = True
    for r in range(1, a.rmax + 1):
        t = Timer(); summary, lines = run([exe, "S", r]); summary["seconds"] = t()
        good = summary["s1_failures"] == 0 and summary["leaves"] == 0 and not lines
        if r in EXPECTED_S1:
            good &= all(summary[x] == y for x, y in EXPECTED_S1[r].items())
        ok &= good; res["s1"][r] = summary
        print(f"   r={r}: {summary['nodes']} nodes, {summary['feasibility_tests']} exact feasibility tests, "
              f"{summary['infeasibility_certificates']} verified infeasibility certificates, failures {summary['s1_failures']} ({summary['seconds']} s)", flush=True)
    checks["s1"] = report(f"(S1) for every sequence with r <= {a.rmax} (attaining structures, exact certificates)", ok)

    if a.r10_chunks:
        t = Timer(); C = a.r10_chunks
        idx = [int(x) for x in a.r10_only.split(",")] if a.r10_only else list(range(C))
        jobs = [[exe, "S", 10, 6, i, C] for i in idx]
        with Pool(max(1, min(a.workers, len(jobs)))) as pool:
            outs = pool.map(run, jobs)
        ok = all(s["s1_failures"] == 0 and s["leaves"] == 0 and not l for s, l in outs)
        if C == 6:
            ok &= all((s["nodes"], s["feasibility_tests"], s["infeasibility_certificates"], s["connectivity_pruned"]) == EXPECTED_R10[i]
                      for i, (s, _) in zip(idx, outs))
        res["s1"]["10"] = [s for s, _ in outs]
        checks["s1_r10"] = report(f"(S1) for r = 10, chunks {idx} of {C} (split at level 6)", ok, f"({t()} s)")

    t = Timer(); ok = True
    for r in range(1, a.catmax + 1):
        summary, lines = run([exe, "C", r])
        recs = [AS.parse_c_line(x) for x in lines if x.startswith("AS ")]
        good = (summary["leaves"], summary["several_classes"], summary["all_degrees_at_least_two"]) == EXPECTED_CAT[r] and len(recs) == summary["leaves"]
        good &= not any(x.startswith("S1 FAIL") for x in lines)
        if a.write_data:
            os.makedirs(os.path.dirname(cat_path(r)), exist_ok=True)
            write_jsonl(cat_path(r), [{"levels": [[list(p) for p in lev] for lev in s], "witness": w} for s, w in recs], compress=r >= 6)
        dep = load_catalogue(r)
        chk = AS.check_catalogue(dep, r)
        good &= {tuple(tuple(tuple(p) for p in lev) for lev in s) for s, _ in dep} == {s for s, _ in recs}
        good &= chk["records"] == chk["distinct"] == EXPECTED_CAT[r][0] and chk["witness_failures"] == 0 and chk["s1_failures"] == 0
        good &= (chk["several_classes"], chk["all_degrees_at_least_two"]) == EXPECTED_CAT[r][1:]
        if r in BRUTE:
            bf = AS.brute_force_structures(r, BRUTE[r]); cat = {tuple(tuple(tuple(p) for p in lev) for lev in s) for s, _ in dep}
            chk["brute_force_structures"] = len(bf); chk["brute_force_missing"] = len(bf - cat)
            good &= not (bf - cat)
        ok &= good; res["catalogue"][r] = {"summary": summary, "independent_check": chk}
        print(f"   catalogue r={r}: {summary['leaves']} structures, {summary['several_classes']} with several classes; deposited catalogue rechecked: {chk}", flush=True)
    checks["catalogue"] = report(f"catalogue of attaining structures r <= {a.catmax}: counts, deposited data, independent witness check, (S1), brute-force completeness",
                                 ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "attaining_structures")


if __name__ == "__main__":
    main()
