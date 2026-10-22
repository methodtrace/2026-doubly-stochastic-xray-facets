#!/usr/bin/env python3
"""Statement (S1) for every real sequence with r <= 9 steps: exact verification of the certificate trees.

Claim (computer-assisted, thm:S1cert): for every real sequence E_0, ..., E_r with 1 <= r <= 9 (and r = 10, see
check_s1_depth10.py), if the boundary co-level classes number at least two, some q < r has at most one partner.
In particular (S1) holds for every composition with at most 10 steps, for every k. (thm:S1 uses (S1) for g and
g*, so this gives prob:B only for k <= 21, inside the range already covered by cor:probB; it extends the
exhaustive (S1) evidence of sec:census from integer compositions to all real sequences.)

Method: see src/dsxray/s1_certificates.py. Each tree partitions the mixed-radix tree of attaining-set choices into
terminals of four kinds (connectivity cut, solitary-index cut, vanishing strict row, Gordan multiplier); each is
checked in exact integer arithmetic, with no use of the program that produced the tree.

Data: data/s1/s1_certificates_rN.txt.gz (r = 1..8), data/s1/s1_certificates_r9.txt.xz, and
data/s1/s1_certificates.json (terminal counts, SHA-256 and size of each uncompressed tree).
Terminal records: 1, 1, 7, 75, 893, 9,293, 104,157, 1,076,939, 11,707,545 for r = 1..9 (12,898,911 in all).

Checks:
  1. every record of each tree for r <= --rmax: interval partition of the whole tree, each certificate, summary
     counts, no counterexample leaf; counts and SHA-256 equal to data/s1/s1_certificates.json;
  2. end-to-end test of the encoding: for --probes random integer sequences per r (small values, so that ties and
     several classes are frequent) the terminal containing the actual attaining structure is a connectivity or
     solitary-index cut whose claim holds for the sequence (attaining_structures.n_classes, .partners), and (S1)
     holds for it (attaining_structures.s1_holds);
  3. (with a C compiler, unless --no-generator) the generator tools/c/s1_certificate_tree.c reproduces the trees for
     r <= --rmax and the r = 10 domain 1_1_3_3_7_7 byte for byte (uncompressed) (SHA-256 of data/s1/s1_certificates.json
     and data/s1/s1_depth10_trees.sha256).
Default: r <= 7 (about 5 s). --full: r <= 9 (about 7 minutes on one core, dominated by r = 9).
"""
import argparse
import os
import random

import hashlib
import shutil
import subprocess

from common import DATA, Timer, build_c, finish, out_dir, read_json, report

from dsxray import attaining_structures as AS
from dsxray import s1_certificates as S1


def tree_path(r):
    return os.path.join(DATA, "s1", f"s1_certificates_r{r}.txt." + ("xz" if r == 9 else "gz"))


def probes(r, n, rng):
    out = []
    for i in range(n):
        lim = 2 + i % 5
        E = [0, 0] + [rng.randint(-lim, lim) for _ in range(r - 1)]
        if i % 2:
            E = [x + 3 * p for p, x in enumerate(E)]
        out.append(E)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/s1_certificates")
    ap.add_argument("--rmax", type=int, default=7)
    ap.add_argument("--probes", type=int, default=3000)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--no-generator", action="store_true",
                    help="do not compile and rerun the C generator (check the deposited certificates only)")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 9)
    out = out_dir(a.out); res = {}; checks = {}
    meta = read_json(os.path.join(DATA, "s1", "s1_certificates.json"))
    rng = random.Random(29); total = 0
    for r in range(1, a.rmax + 1):
        t = Timer(); P = probes(r, a.probes, rng)
        s1ok = all(AS.s1_holds(AS.structure_of(E), r) for E in P)
        multi = sum(AS.n_classes(AS.structure_of(E), r) >= 2 for E in P)
        try:
            s = S1.verify_tree(tree_path(r), r, (), P, progress=2000000 if r == 9 else None); err = None
        except AssertionError as e:
            s, err = None, str(e)
        m = meta[str(r)]
        ok = (s is not None and s1ok and s["probes_checked"] == len(P)
              and (s["terminal_records"], s["certificate_types"], s["uncompressed_sha256"], s["uncompressed_bytes"])
              == (m["terminal_records"], m["certificate_types"], m["uncompressed_sha256"], m["uncompressed_bytes"]))
        res[f"r{r}"] = dict(s, probes_with_several_classes=multi) if s else {"error": err}
        total += s["terminal_records"] if s else 0
        checks[f"tree_r{r}"] = report(f"(S1) certificate tree, r = {r}: complete and every terminal certified", ok,
                                      f"{s['terminal_records'] if s else err} terminals, {len(P)} probes ({multi} with several classes) ({t()} s)")
    res["terminal_records_total"] = total
    if not a.no_generator and (shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")):
        t = Timer(); exe = build_c("s1_certificate_tree", os.path.join(out, "bin")); ok = True
        for r in range(1, a.rmax + 1):
            h = hashlib.sha256()
            with subprocess.Popen([exe, "S", str(r), "-1", "0", "1", "1", "1", "0", "_"], stdout=subprocess.PIPE) as proc:
                for chunk in iter(lambda: proc.stdout.read(1 << 22), b""):
                    h.update(chunk)
            ok &= proc.returncode == 0 and h.hexdigest() == meta[str(r)]["uncompressed_sha256"]
        out_bytes = subprocess.run([exe, "S", "10", "-1", "0", "1", "1", "1", "0", "1,1,3,3,7,7"], capture_output=True, check=True).stdout
        with open(os.path.join(DATA, "s1", "s1_depth10_trees.sha256")) as f:
            want = dict(reversed(line.split()) for line in f if line.strip())["1_1_3_3_7_7.txt"]
        ok &= hashlib.sha256(out_bytes).hexdigest() == want
        checks["generator"] = report("the C generator reproduces the trees byte for byte", ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "s1_certificates")


if __name__ == "__main__":
    main()
