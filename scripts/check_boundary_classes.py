#!/usr/bin/env python3
"""Statement (S1), lem:R and cor:probB (total support implies ray for staircases with k <= 30).

  1. (S1): whenever the boundary steps form at least two classes, some q < r is solitary (at most one attaining
     partner). Python check for k <= --kmax; with --c-kmax K the C program tools/c/boundary_classes.c checks every
     composition with 4 <= k <= K and 1 <= r <= k-2. The paper's range k <= 30 has 1,073,741,762 compositions, of which
     483,464,468 have at least two classes (cor:probB, tab:evidence); these totals are checked on the deposited counts
     (data/summary/boundary_classes.json) and, with --c-kmax 30, recomputed. --lifts adds (S1) on pattern lifts.
  2. lem:R: every composition with a solitary index lacks total support, from the tight graph computed directly
     from v = u box u (k <= --lemma-r-kmax).
  3. The composition (1,1,4,4,7,4,4,1,6,5) (k = 37) of sec:s_census-proofs and two relatives have a tight
     permutation and no total support, two classes, dim Z = 3, and satisfy (S1); the k = 22 example
     (1,1,5,3,2,2,3,5) of sec:SM-square has one class, total support and is a ray (see also check_examples.py).
  4. Additional, not used in the paper: statement (H+) (a dichotomy for compositions with
     several classes) holds for k <= --kmax (Python) and k <= --c-hplus-kmax (C, tools/c/hplus_statement.c) but fails
     at k = 37 (the compositions of item 3).
Defaults: Python k <= 16, lem:R k <= 12, no C. --full: C for (S1) k <= 30 (about 40 minutes), (H+) k <= 31,
pattern lifts, lem:R k <= 14 (about 1.5 hours on one core in total).
"""
import argparse
import json
import os
import shutil
import subprocess
import tempfile

from common import DATA, ROOT, Timer, finish, out_dir, read_json, report

from dsxray import boundary_classes as BC
from dsxray.level_graph import hplus_counts, several_classes_without_deficient_class
from dsxray.staircases import staircases

COUNTEREXAMPLES = [(1, 1, 4, 4, 7, 4, 4, 1, 6, 5), (2, 1, 4, 4, 7, 4, 4, 1, 6, 5), (1, 7, 10, 10, 13, 10, 10, 7, 12, 11)]
ONE_CLASS_EXAMPLE = (1, 1, 5, 3, 2, 2, 3, 5)
LIFTS = [(2, 13), (3, 12), (4, 10), (5, 8), (7, 7)]


def build(src, workdir):
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if not cc:
        raise SystemExit("a C compiler (cc) is needed for --c-kmax / --c-hplus-kmax / --lifts")
    exe = os.path.join(workdir, os.path.basename(src)[:-2])
    subprocess.run([cc, "-O2", "-o", exe, src], check=True)
    return exe


def run_c(exe, args):
    p = subprocess.run([exe] + [str(x) for x in args], check=True, capture_output=True, text=True)
    return [json.loads(line) for line in p.stdout.splitlines() if line.startswith("{")], \
           [line for line in p.stdout.splitlines() if not line.startswith("{")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/boundary_classes")
    ap.add_argument("--kmax", type=int, default=16)
    ap.add_argument("--lemma-r-kmax", type=int, default=12)
    ap.add_argument("--c-kmax", type=int, default=0)
    ap.add_argument("--c-hplus-kmax", type=int, default=0)
    ap.add_argument("--lifts", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.c_kmax = max(a.c_kmax, 30); a.c_hplus_kmax = max(a.c_hplus_kmax, 31); a.lifts = True
        a.lemma_r_kmax = max(a.lemma_r_kmax, 14)
    dep = read_json(os.path.join(DATA, "summary", "boundary_classes.json"))
    dep_s1 = {r["k"]: r for r in dep["s1_compositions"]}
    dep_h = {r["k"]: r for r in dep["hplus_compositions"]}
    res = {}; checks = {}

    # 1 and 4: Python
    t = Timer(); s1 = [BC.s1_counts(k) for k in range(4, a.kmax + 1)]
    res["s1_python"] = s1
    checks["s1_python"] = report(f"(S1) on all compositions, k <= {a.kmax} (Python)",
                                 all(r["s1_failures"] == 0 and r == dep_s1[r["k"]] for r in s1), f"({t()} s)")
    t = Timer(); hp = [hplus_counts(k) for k in range(4, a.kmax + 1)]
    res["hplus_python"] = hp
    checks["hplus_python"] = report(f"(H+) on all compositions, k <= {a.kmax} (Python)",
                                     all("FAIL" not in r["counts"] and sum(r["counts"].values()) == dep_h[r["k"]]["compositions"]
                                         and r["counts"].get("connected", 0) == dep_h[r["k"]]["connected"] for r in hp), f"({t()} s)")

    if a.kmax >= 16:
        t = Timer(); sev = nodef = 0
        for k in range(4, 17):
            for g in staircases(k):
                s, nd = several_classes_without_deficient_class(g); sev += s; nodef += nd
        res["several_classes_k16"] = {"several_classes": sev, "no_deficient_class": nodef}
        checks["balanced_case_k16"] = report("k <= 16: 9,299 compositions with several classes, 4,842 without a deficient class",
                                             (sev, nodef) == (9299, 4842), f"({t()} s)")

    tot = sum(r["compositions"] for r in dep["s1_compositions"] if 4 <= r["k"] <= 30)
    sev = sum(r["several_classes"] for r in dep["s1_compositions"] if 4 <= r["k"] <= 30)
    fl = sum(r["s1_failures"] for r in dep["s1_compositions"] if 4 <= r["k"] <= 30)
    res["s1_deposited_totals_k4_to_30"] = {"compositions": tot, "several_classes": sev, "failures": fl}
    checks["s1_deposited_totals"] = report("deposited (S1) counts for 4 <= k <= 30: 1,073,741,762 compositions, "
                                           "483,464,468 with several classes, no failure",
                                           (tot, sev, fl) == (1073741762, 483464468, 0) and tot == sum(2 ** (k - 1) - 2 for k in range(4, 31)))

    # 2: lem:R
    t = Timer(); lr = BC.lemma_r_counts(4, a.lemma_r_kmax)
    res["lemma_r"] = dict(lr, k_range=[4, a.lemma_r_kmax])
    ok = lr["of_which_total_support"] == 0
    if a.lemma_r_kmax == 14:
        ok = ok and lr["with_solitary_index"] == 5272
    checks["lemma_r"] = report(f"lem:R, k <= {a.lemma_r_kmax}", ok, f"{lr} ({t()} s)")

    # 3: (H+) counterexamples
    desc = [BC.describe(g) for g in COUNTEREXAMPLES]
    res["hplus_counterexamples"] = desc
    checks["hplus_false"] = report("(H+) fails at k = 37, 38, 91 (no total support, two classes, rays, (S1) holds)",
                                   all(d["hplus"] == "FAIL" and not d["total_support"] and d["unused_tight_edges"]
                                       and len(d["classes"]) == 2 and d["dim_z"] == 3 and d["s1_holds"] for d in desc),
                                   f"k=37: {desc[0]['unused_tight_edges']} of {desc[0]['tight_edges']} tight edges in no tight permutation")
    d = BC.describe(ONE_CLASS_EXAMPLE); res["one_class_example"] = d
    checks["k22_example"] = report("k = 22 example: one class, total support, ray",
                                   len(d["classes"]) == 1 and d["total_support"] and d["dim_z"] == 3)

    # C programs
    with tempfile.TemporaryDirectory() as tmp:
        if a.c_kmax or a.lifts:
            exe = build(os.path.join(ROOT, "tools", "c", "boundary_classes.c"), tmp)
        if a.c_kmax:
            t = Timer(); rows, fails = run_c(exe, ["A", 4, a.c_kmax]); res["s1_c"] = rows
            checks["s1_c"] = report(f"(S1) on all compositions, k <= {a.c_kmax} (C)",
                                    not fails and all(r["s1_failures"] == 0 and r == dep_s1[r["k"]] for r in rows),
                                    f"{sum(r['compositions'] for r in rows)} compositions ({t()} s)")
        if a.lifts:
            t = Timer(); allrows = []
            for V, R in LIFTS:
                rows, fails = run_c(exe, ["P", V, R]); allrows += rows
                checks[f"lifts_V{V}"] = report(f"(S1) on pattern lifts V={V}, r <= {R}", not fails and all(r["s1_failures"] == 0 for r in rows))
            res["s1_lifts"] = allrows
        if a.c_hplus_kmax:
            exe2 = build(os.path.join(ROOT, "tools", "c", "hplus_statement.c"), tmp)
            t = Timer(); rows, fails = run_c(exe2, [4, a.c_hplus_kmax]); res["hplus_c"] = rows
            checks["hplus_c"] = report(f"(H+) on all compositions, k <= {a.c_hplus_kmax} (C)",
                                       not fails and all(r["FAIL"] == 0 and r == dep_h[r["k"]] for r in rows), f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "boundary_classes")


if __name__ == "__main__":
    main()
