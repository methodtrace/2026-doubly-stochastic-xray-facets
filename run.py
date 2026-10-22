#!/usr/bin/env python3
"""Run the test suite and every reproduction script at one of four tiers.

  python3 run.py                     quick tier: tests and every script at its small default range (about 7 to 25
                                     minutes on 2 cores, depending on the machine; needs pytest and a C compiler)
  python3 run.py --tier certificates verify the shipped certificates and data without recomputing them: the
                                     checksums of MANIFEST.sha256, the certified classification of the deposited
                                     rays, the E-flat sign-cell streams r <= 7, the (S1) trees r <= 9, the r = 10
                                     index (and the r = 10 domain trees with --deposit)
                                     (about 11 minutes on 2 cores; no C compiler and no pytest needed)
  python3 run.py --tier medium       larger ranges (a few hours on 2 cores)
  python3 run.py --tier full         the ranges of the full tier (see README.md for the cost of each step; about
                                     150 core-hours, dominated by the staircase sweep at k = 23, 24)
  python3 run.py --tier paper        every range stated in the paper, including the walks at k = 12, the line
                                     enumeration at k = 10 and the regeneration of every certificate tree, as
                                     independent pieces (--list prints them, --piece ID runs one; a rerun skips the
                                     pieces that passed). OUT/pieces/<id>/ holds each piece's results, log and timing.

Options: --workers N (default 2), --out DIR (default results), --only NAME[,NAME...] (run selected steps),
--skip-tests, --deposit PATH (the separate file s1_depth10_certificates.zip, or its unpacked directory; passed to the
s1_depth10 step). Each step writes OUT/<step>/<name>.json; OUT/run_summary.json lists every step with its status and,
where a step checks less than the paper's range at this tier, its scope (for example "index only").
Exit status 0 iff every step passed.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from math import comb

ROOT = os.path.dirname(os.path.abspath(__file__))

STEPS = ["tests", "facet_counts", "census", "shapes", "vertex_walk", "staircases", "nonincreasing_tails", "three_step",
         "hall_theorem", "boundary_classes", "cells", "cells_depth7", "staircase_subclasses", "attaining_structures",
         "s1_certificates", "s1_depth10", "facet_criterion", "facets_without_total_support", "search_total_support",
         "facet_fibers", "padding", "refinement", "refinement_descent", "two_point_lifts", "families",
         "family_distinctness", "unbounded_height", "quadratic_family", "height_two_rays", "binary_increments",
         "examples"]

SCRIPTS = {
    "facet_counts": "reproduce_facet_counts.py", "census": "reproduce_census.py", "shapes": "reproduce_shapes.py",
    "vertex_walk": "crosscheck_vertex_walk.py", "staircases": "check_staircase_conjectures.py",
    "three_step": "check_three_step_theorem.py", "hall_theorem": "check_hall_theorem.py",
    "boundary_classes": "check_boundary_classes.py", "cells": "check_cells.py", "padding": "check_padding.py",
    "families": "check_families.py", "nonincreasing_tails": "check_nonincreasing_tails.py",
    "family_distinctness": "check_family_distinctness.py", "examples": "check_examples.py",
    "attaining_structures": "enumerate_attaining_structures.py",
    "facet_criterion": "check_facet_criterion.py", "facets_without_total_support": "check_facets_without_total_support.py",
    "search_total_support": "search_total_support.py", "refinement": "check_refinement.py",
    "unbounded_height": "check_unbounded_height.py",
    "cells_depth7": "check_cells_depth7.py", "staircase_subclasses": "check_staircase_subclasses.py",
    "s1_certificates": "check_s1_certificates.py", "s1_depth10": "check_s1_depth10.py",
    "s1_regenerate": "regenerate_s1_depth10.py",
    "facet_fibers": "check_facet_fibers.py", "quadratic_family": "check_quadratic_family.py",
    "refinement_descent": "check_refinement_descent.py", "two_point_lifts": "check_two_point_lifts.py",
    "height_two_rays": "check_height_two_rays.py", "binary_increments": "check_binary_increments.py"}

QUICK = {"facets_without_total_support": ["--search"]}

MEDIUM = {
    "facet_counts": ["--kmax", "12", "--walk-kmax", "10", "--line-kmax", "9"],
    "vertex_walk": ["--kmax", "10", "--cdd-kmax", "10", "--dd-kmax", "8"],
    "staircases": ["--kmax", "18"],
    "nonincreasing_tails": ["--full"], "family_distinctness": ["--full"],
    "three_step": ["--proof-kmax", "20", "--kmax", "30"],
    "hall_theorem": ["--kmax", "16", "--brute-kmax", "11", "--proof-kmax", "12", "--sizes-kmax", "13"],
    "boundary_classes": ["--kmax", "16", "--lemma-r-kmax", "14", "--c-kmax", "26", "--c-hplus-kmax", "26", "--lifts"],
    "cells": ["--rmax", "6", "--nmax", "20", "--kmax", "14"],
    "padding": ["--full"], "families": ["--full"],
    "attaining_structures": ["--rmax", "8"], "facet_criterion": ["--samples", "100"],
    "facets_without_total_support": ["--search", "--search-kmax", "15", "--rmax", "12"],
    "search_total_support": ["--inc-kmax", "15"], "refinement": ["--full"], "unbounded_height": ["--bmax", "25"],
    "staircase_subclasses": ["--full"], "s1_certificates": ["--rmax", "8"],
    "facet_fibers": ["--full"], "quadratic_family": ["--bmax", "10", "--classify-bmax", "6", "--kmax", "24"],
    "refinement_descent": ["--full"], "two_point_lifts": ["--full"], "height_two_rays": ["--kmax", "30"],
    "binary_increments": ["--full"]}

# quick-tier steps that compile a C program unconditionally (the other steps skip their C parts without a compiler)
CC_STEPS = {"attaining_structures", "search_total_support", "facets_without_total_support"}

# certificate tier: (step, extra arguments); every step reads deposited certificates or data and recomputes nothing
# that the certificates already decide (no searches, no enumerations, no C programs)
CERTIFICATES = [("manifest", None),
                ("facet_counts", ["--walk-kmax", "4", "--profiles-kmax", "1", "--line-kmax", "8"]),
                ("cells_depth7", ["--full", "--no-generator"]),
                ("s1_certificates", ["--full", "--no-generator"]),
                ("s1_depth10", [])]

WORKER_STEPS = {"facet_counts", "vertex_walk", "staircases", "nonincreasing_tails", "three_step", "hall_theorem",
                "attaining_structures", "facets_without_total_support", "search_total_support", "s1_depth10",
                "s1_regenerate"}

DEPOSIT_ZIP = "s1_depth10_certificates.zip"


def _balanced(weights, target):
    """Split range(len(weights)) into consecutive runs whose weights stay near target."""
    runs, start, acc = [], 0, 0
    for i, w in enumerate(weights):
        if acc and acc + w > target:
            runs.append((start, i)); start, acc = i, 0
        acc += w
    runs.append((start, len(weights)))
    return runs


def paper_pieces():
    """The paper tier as independent pieces: (id, step, arguments, dependencies). Every range stated in the paper is
    reached; pieces can run in any order that respects the dependencies, on one machine or several, and the result
    of a piece does not depend on the machine or on the number of workers."""
    P = []

    def add(pid, step, args, deps=()):
        P.append({"id": pid, "step": step, "args": list(args), "deps": list(deps)})

    add("tests", "tests", [])
    add("facet-counts", "facet_counts", ["--full", "--write-facets"])
    add("facet-walk-k12", "facet_counts", ["--kmax", "4", "--walk-kmin", "12", "--walk-kmax", "12",
                                           "--profiles-kmax", "1", "--line-kmax", "8"])
    add("line-enumeration-k10", "facet_counts", ["--kmax", "4", "--walk-kmax", "4", "--profiles-kmax", "1",
                                                 "--line-kmin", "10", "--line-kmax", "10"])
    add("census", "census", ["--full"])
    add("shapes", "shapes", ["--full"])
    add("vertex-walk", "vertex_walk", ["--full", "--require-cdd"])
    add("vertex-walk-exact-k12", "vertex_walk", ["--kmin", "12", "--kmax", "12", "--cdd-kmax", "4", "--dd-kmax", "3"])
    add("vertex-walk-cdd-k11-k12", "vertex_walk", ["--kmax", "4", "--cdd-kmin", "11", "--cdd-kmax", "12", "--dd-kmax", "3",
                                                   "--require-cdd"])
    stair = ["staircases-k3-k12"]
    add("staircases-k3-k12", "staircases", ["--kmax", "12"])
    for k in range(13, 25):
        sizes = [comb(k - 1, r) for r in range(1, k - 1)]
        for a, b in _balanced(sizes, 1_500_000):
            pid = f"staircases-k{k}" if (a, b) == (0, len(sizes)) else f"staircases-k{k}-r{a + 1}-r{b}"
            add(pid, "staircases", ["--k", str(k), "--r", f"{a + 1}:{b}", "--no-small"]); stair.append(pid)
    add("staircases-merge", "staircases", ["--merge"] + [f"@{p}" for p in stair], deps=stair)
    add("nonincreasing-tails", "nonincreasing_tails", ["--full"])
    add("three-step", "three_step", ["--full"])
    add("hall-theorem", "hall_theorem", ["--full"])
    add("boundary-classes", "boundary_classes", ["--full"])
    add("cells", "cells", ["--full"])
    add("cells-depth7", "cells_depth7", ["--full"])
    add("staircase-subclasses", "staircase_subclasses", ["--full"])
    add("attaining-structures", "attaining_structures", ["--full"])
    for name, chunks in (("a", "0,1"), ("b", "2,3"), ("c", "4,5")):
        add(f"attaining-structures-r10-{name}", "attaining_structures",
            ["--rmax", "1", "--catmax", "1", "--r10-chunks", "6", "--r10-only", chunks])
    add("s1-certificates", "s1_certificates", ["--full"])
    with open(os.path.join(ROOT, "data", "s1", "s1_depth10_domains.json")) as f:
        sizes = [d["uncompressed_bytes"] for d in json.load(f)["domain_list"]]
    regen = []
    for n, (a, b) in enumerate(_balanced(sizes, 600_000_000), 1):
        pid = f"s1-regenerate-{n:02d}"; regen.append(pid)
        add(pid, "s1_regenerate", ["--domains", f"{a}:{b}"])
    add("s1-deposit", "s1_depth10", ["--make-deposit", ",".join(f"@{p}/trees" for p in regen), f"@s1-deposit/{DEPOSIT_ZIP}",
                                     "--no-verify"], deps=regen)
    for pid, step in (("facet-criterion", "facet_criterion"), ("facets-without-total-support", "facets_without_total_support"),
                      ("search-total-support", "search_total_support"),
                      ("padding", "padding"), ("refinement", "refinement"), ("refinement-descent", "refinement_descent"),
                      ("two-point-lifts", "two_point_lifts"), ("families", "families"),
                      ("family-distinctness", "family_distinctness"), ("unbounded-height", "unbounded_height"),
                      ("quadratic-family", "quadratic_family"), ("height-two-rays", "height_two_rays"),
                      ("binary-increments", "binary_increments"), ("examples", "examples")):
        add(pid, step, ["--full"])
    add("facet-fibers", "facet_fibers", ["--full", "--mmax", "6"])
    return P


def _resolve(arg, out):
    """@ID -> OUT/pieces/ID (a piece's folder); @@NAME -> OUT/NAME; inside comma lists as well."""
    def one(x):
        if x.startswith("@@"):
            return os.path.join(out, x[2:])
        if x.startswith("@"):
            return os.path.join(out, "pieces", x[1:])
        return x
    return ",".join(one(x) for x in arg.split(",")) if "@" in arg else arg


def _env():
    return dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src") + os.pathsep + os.environ.get("PYTHONPATH", ""),
                PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")


def _peak_bytes(ru):
    return ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024


def run_piece(piece, out, workers):
    """Run one piece into OUT/pieces/<id>; record status, wall and CPU time and peak memory in piece.json there."""
    import platform
    import resource
    pdir = os.path.join(out, "pieces", piece["id"])
    for dep in piece["deps"]:
        try:
            with open(os.path.join(out, "pieces", dep, "piece.json")) as f:
                if json.load(f)["status"] != "pass":
                    raise OSError
        except (OSError, KeyError, ValueError):
            raise SystemExit(f"piece {piece['id']} needs piece {dep}, which has not passed yet")
    os.makedirs(pdir, exist_ok=True)
    if piece["step"] == "tests":
        cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", os.path.join(ROOT, "tests")]
    else:
        cmd = [sys.executable, os.path.join(ROOT, "scripts", SCRIPTS[piece["step"]]), "--out", pdir]
        if piece["step"] in WORKER_STEPS:
            cmd += ["--workers", str(workers)]
        cmd += [_resolve(x, out) for x in piece["args"]]
    print(f"\n=== piece {piece['id']}: {' '.join(os.path.relpath(c, ROOT) if c.startswith(ROOT) else c for c in cmd[1:])}",
          flush=True)
    t = time.time(); r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    with open(os.path.join(pdir, "log.txt"), "w") as log:
        proc = subprocess.Popen(cmd, cwd=ROOT, env=_env(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            sys.stdout.write(line); log.write(line)
        rc = proc.wait()
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    rec = {"id": piece["id"], "step": piece["step"], "args": piece["args"], "status": "pass" if rc == 0 else "FAIL",
           "wall_seconds": float(f"{time.time() - t:.1f}"),
           "cpu_seconds": float(f"{(r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime):.1f}"),
           "peak_memory_bytes_largest_process": _peak_bytes(r1), "workers": workers,
           "machine": {"system": platform.system(), "machine": platform.machine(), "processor": platform.processor(),
                       "cpus": os.cpu_count(), "python": platform.python_version()},
           "finished": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    with open(os.path.join(pdir, "piece.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(f"=== piece {piece['id']}: {rec['status']} ({rec['wall_seconds']} s wall, {rec['cpu_seconds']} s CPU)", flush=True)
    return rc == 0


def paper_main(a):
    pieces = paper_pieces()
    if a.list:
        print(json.dumps(pieces, indent=1))
        return 0
    byid = {p["id"]: p for p in pieces}
    if a.piece:
        bad = [x for x in a.piece.split(",") if x not in byid]
        if bad:
            raise SystemExit(f"unknown piece(s): {', '.join(bad)} (python3 run.py --tier paper --list)")
        return 0 if all([run_piece(byid[x], a.out, a.workers) for x in a.piece.split(",")]) else 1
    ok = True
    for p in pieces:
        done = os.path.join(a.out, "pieces", p["id"], "piece.json")
        if os.path.exists(done):
            with open(done) as f:
                if json.load(f)["status"] == "pass":
                    continue
        ok &= run_piece(p, a.out, a.workers)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["quick", "certificates", "medium", "full", "paper"], default="quick")
    ap.add_argument("--list", action="store_true", help="paper tier: print the pieces as JSON and exit")
    ap.add_argument("--piece", default="", help="paper tier: run these pieces (comma-separated ids) and exit")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(ROOT, "results"))
    ap.add_argument("--only", default="")
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--deposit", default="", help="s1_depth10_certificates.zip (or its directory) for the r = 10 trees")
    a = ap.parse_args()
    if sys.version_info[:2] < (3, 12):
        raise SystemExit("Python 3.12 or newer is required")
    if a.tier == "paper" or a.piece or a.list:
        return paper_main(a)
    all_steps = [s for s, _ in CERTIFICATES] if a.tier == "certificates" else STEPS
    steps = [s for s in all_steps if not a.only or s in a.only.split(",")]
    if a.skip_tests and "tests" in steps:
        steps.remove("tests")
    missing = []
    if "tests" in steps and importlib.util.find_spec("pytest") is None:
        missing.append("pytest (python3 -m pip install pytest), or pass --skip-tests")
    needs_cc = a.tier in ("medium", "full") or (a.tier == "quick" and CC_STEPS & set(steps))
    if needs_cc and not (shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")):
        missing.append("a C compiler (cc, gcc or clang); without one use --tier certificates")
    if missing:
        raise SystemExit("missing prerequisites for this tier:\n  " + "\n  ".join(missing))
    os.makedirs(a.out, exist_ok=True)
    env = _env()
    summary = []; t_all = time.time()
    cert_args = dict(CERTIFICATES)
    for s in steps:
        scope = ""
        if s == "manifest":
            print("\n=== manifest: SHA-256 of every file listed in MANIFEST.sha256", flush=True)
            t = time.time(); bad = []; n = 0
            with open(os.path.join(ROOT, "MANIFEST.sha256")) as f:
                for line in f:
                    if not line.strip():
                        continue
                    want, name = line.rstrip("\n").split("  ", 1); n += 1
                    try:
                        with open(os.path.join(ROOT, name), "rb") as g:
                            got = hashlib.file_digest(g, "sha256").hexdigest()
                    except OSError:
                        got = None
                    if got != want:
                        bad.append(name)
            for name in bad:
                print(f"   checksum mismatch or missing: {name}")
            print(f"   {n - len(bad)} of {n} files OK", flush=True)
            summary.append({"step": s, "status": "FAIL" if bad else "pass", "seconds": float(f"{time.time() - t:.1f}")})
            print(f"=== {s}: {summary[-1]['status']} ({summary[-1]['seconds']} s)", flush=True)
            continue
        if s == "tests":
            cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", os.path.join(ROOT, "tests")]
        else:
            cmd = [sys.executable, os.path.join(ROOT, "scripts", SCRIPTS[s]), "--out", os.path.join(a.out, s)]
            if s in WORKER_STEPS:
                cmd += ["--workers", str(a.workers)]
            if a.tier == "certificates":
                cmd += cert_args[s]
            elif a.tier == "quick":
                cmd += QUICK.get(s, [])
            elif a.tier == "medium":
                cmd += MEDIUM.get(s, [])
            elif a.tier == "full":
                cmd += ["--full"]
            if s == "s1_depth10":
                if a.deposit:
                    cmd += ["--deposit", a.deposit]
                else:
                    scope = "index only (r = 10 domain trees not checked: pass --deposit)"
        print(f"\n=== {s}: {' '.join(os.path.relpath(c, ROOT) if c.startswith(ROOT) else c for c in cmd[1:])}", flush=True)
        t = time.time(); rc = subprocess.call(cmd, cwd=ROOT, env=env)
        summary.append({"step": s, "status": "pass" if rc == 0 else "FAIL", "seconds": float(f"{time.time() - t:.1f}")})
        if scope:
            summary[-1]["scope"] = scope
        print(f"=== {s}: {'pass' if rc == 0 else 'FAIL'} ({summary[-1]['seconds']} s){'; ' + scope if scope else ''}", flush=True)
    ok = all(x["status"] == "pass" for x in summary)
    with open(os.path.join(a.out, "run_summary.json"), "w") as f:
        json.dump({"tier": a.tier, "steps": summary, "all_passed": ok}, f, indent=1)
    print(f"\n{'ALL STEPS PASSED' if ok else 'SOME STEPS FAILED'} ({a.tier} tier, {time.time() - t_all:.0f} s)")
    for x in summary:
        print(f"  {x['status']:4s}  {x['seconds']:8.1f} s  {x['step']}" + (f"  [{x['scope']}]" if "scope" in x else ""))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
