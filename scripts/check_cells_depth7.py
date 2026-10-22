#!/usr/bin/env python3
"""Statement E-flat for r <= 7 steps: exact verification of the complete sign-cell certificates.

Claim (computer-assisted): statement E-flat holds for every tail (g_1, ..., g_r) of positive integers with r <= 7.
This is thm:Eflat6; with thm:Eflat it gives cor:probA: prob:A holds for every staircase with r <= 7 or r >= k - 8,
and the rigid-corner statement conj:E holds for every staircase with r <= 7. Table of nonempty cells (r, forms,
cells, Hall fails, K* empty, E-flat holds, fails):
  r = 1..6 as in tab:cells (2, 6, 26, 198, 2,298, 45,014 cells), and r = 7: 34 forms, 1,276,978 cells,
  976,412 / 16,682 / 283,884 / 0.
Method (see src/dsxray/sign_certificates.py): by the two-cell reduction, each realizable full sign vector sigma of
the comparison forms gives exactly the two cells C(empty, sigma) and C(J(sigma), sigma); the stream lists the
realizable sign vectors with integer witnesses and every rejected sign prefix with an exact infeasibility
certificate (Gordan multipliers or a vanishing strict row), and the terminal prefixes partition the base-3 sign
tree. No floating point and no Fourier-Motzkin elimination are used; this is independent of check_cells.py.

Data: data/cells/sign_certificates_rN.txt.gz (r = 1..6) and data/cells/sign_certificates_r7.txt.xz (638,489 sign
vectors and 7,048,460 rejected prefixes, 4,834,874 by Gordan multipliers and 2,213,586 by vanishing rows; 413,497,888 bytes uncompressed, SHA-256 in
data/cells/sign_certificates_r7.json).

Checks:
  1. r = 1..6 (1, 3, 13, 99, 1,149, 22,507 realizable sign vectors): every record of the stream (witnesses, certificates, partition, summary); the expanded cells and their
     verdicts equal data/cells/cells_rN.json exactly; every verdict is also computed by the matching-based routine
     cell_enumeration.eflat_tail (both cells of every sign vector);
  2. r = 6: every integer tail with sum <= --tail-sum lies in a certified cell with the same verdict;
  3. (--full) r = 7: every record of the stream, the counts above, the SHA-256 of the decompressed stream, the
     matching-based verdict on both cells of every realizable sign vector; every integer tail with r = 7 and sum <= 20 (77,520 tails)
     and 20,000 random tails (parts up to 1000, a third of them with forced ties) lie in a certified cell with the
     same verdict as computed directly from the tail;
  4. (with a C compiler, unless --no-generator) the generator tools/c/sign_cell_certificates.c reproduces the
     streams for r <= 6 byte for byte, and with --full also the stream for r = 7 (SHA-256 of the uncompressed output
     equal to data/cells/sign_certificates_r7.json).
Default: r <= 6 (about 10 s). --full: adds r = 7 (about 4 minutes on one core).
"""
import argparse
import hashlib
import lzma
import os
import random
from itertools import combinations

import gzip
import shutil
import subprocess

from common import DATA, Timer, build_c, finish, out_dir, read_json, report

from dsxray import cell_enumeration as CE
from dsxray import sign_certificates as SC

EXPECTED = {1: (2, 0, 2, 0), 2: (6, 0, 6, 0), 3: (26, 0, 20, 6), 4: (198, 12, 90, 96), 5: (2298, 586, 434, 1278),
            6: (45014, 24088, 2654, 18272), 7: (1276978, 976412, 16682, 283884)}
FORMS = {1: 0, 2: 1, 3: 3, 4: 7, 5: 13, 6: 22, 7: 34}
SIGN_VECTORS = {1: 1, 2: 3, 3: 13, 4: 99, 5: 1149, 6: 22507, 7: 638489}


def stream_path(r):
    return os.path.join(DATA, "cells", f"sign_certificates_r{r}.txt." + ("xz" if r == 7 else "gz"))


def tails_with_sum_at_most(r, n):
    """All compositions (g_1..g_r) of positive integers with sum <= n."""
    for tot in range(r, n + 1):
        for cuts in combinations(range(1, tot), r - 1):
            b = (0,) + cuts + (tot,)
            yield tuple(b[i + 1] - b[i] for i in range(r))


def sign_string(t, forms):
    s = []
    for f in forms:
        v = sum(a * b for a, b in zip(f, t))
        s.append("+" if v > 0 else "-" if v < 0 else "0")
    return "".join(s)


def tail_coverage(r, tails, leaves):
    forms = SC.comparison_forms(r); n = 0
    for t in tails:
        sig = sign_string(t, forms)
        v1, J, v2 = leaves[sig]
        units = tuple(i + 1 for i, g in enumerate(t) if g == 1)
        assert units in ((), J), ("tail in an empty cell", t)
        assert SC.verdict_of_tail(t) == (v1 if not units else v2), ("verdict differs", t)
        n += 1
    return n


def cell_counts(res):
    c = res["cell_verdicts"]
    return (res["cells"], c.get("nohall", 0), c.get("ts", 0), c.get("ok", 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/cells_depth7")
    ap.add_argument("--rmax", type=int, default=6)
    ap.add_argument("--tail-sum", type=int, default=16)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--no-generator", action="store_true",
                    help="do not compile and rerun the C generator (check the deposited certificates only)")
    a = ap.parse_args()
    if a.full:
        a.rmax = 7
    out = out_dir(a.out); res = {}; checks = {}
    for r in range(1, min(a.rmax, 6) + 1):
        t = Timer(); leaves = {} if r == 6 else None
        name = f"cells_r{r}.json" if r < 6 else "cells_r6.json.gz"
        dep = read_json(os.path.join(DATA, "cells", name))
        try:
            s = SC.verify_stream(stream_path(r), r, dep, CE.eflat_tail, 1, collect=leaves); err = None
        except AssertionError as e:
            s, err = None, str(e)
        ok = s is not None and s["forms"] == FORMS[r] and cell_counts(s) == EXPECTED[r] and "FAIL" not in s["cell_verdicts"]
        ok = ok and s["sign_vectors"] == SIGN_VECTORS[r]
        res[f"r{r}"] = s or {"error": err}
        checks[f"stream_r{r}"] = report(f"r={r}: complete certificate stream, cells equal to the deposited list", ok,
                                        f"{s['cells'] if s else err} cells ({t()} s)")
        if r == 6 and s is not None:
            t = Timer()
            try:
                n = tail_coverage(6, tails_with_sum_at_most(6, a.tail_sum), leaves); ok = True
            except (AssertionError, KeyError) as e:
                n, ok = str(e), False
            res["r6_tails"] = {"sum_bound": a.tail_sum, "tails": n}
            checks["tails_r6"] = report(f"r=6: every integer tail with sum <= {a.tail_sum} lies in a certified cell with the same verdict",
                                        ok, f"{n} tails ({t()} s)")
    if a.rmax >= 7:
        t = Timer(); meta = read_json(os.path.join(DATA, "cells", "sign_certificates_r7.json"))
        h = hashlib.sha256(); nbytes = 0
        with lzma.open(stream_path(7), "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk); nbytes += len(chunk)
        checks["hash_r7"] = report("r=7: SHA-256 of the decompressed stream", (h.hexdigest(), nbytes) ==
                                   (meta["uncompressed_sha256"], meta["uncompressed_bytes"]), f"{nbytes} bytes ({t()} s)")
        t = Timer(); leaves = {}
        try:
            s = SC.verify_stream(stream_path(7), 7, None, CE.eflat_tail, 1, progress=1000000, collect=leaves); err = None
        except AssertionError as e:
            s, err = None, str(e)
        ok = s is not None and s["forms"] == 34 and cell_counts(s) == EXPECTED[7] and "FAIL" not in s["cell_verdicts"]
        ok = ok and (s["sign_vectors"], s["rejected_prefixes"]) == (meta["sign_vectors"], meta["rejected_prefixes"])
        ok = ok and s["records"] == 7686949 and s["rejection_kinds"] == {"G": 4834874, "Z": 2213586}
        ok = ok and s["second_verdict_leaves"] == s["sign_vectors"]
        res["r7"] = s or {"error": err}
        checks["stream_r7"] = report("r=7: complete certificate stream, 1,276,978 cells, E-flat holds on every cell with K* nonempty",
                                     ok, f"{s['cells'] if s else err} cells ({t()} s)")
        if s is not None:
            t = Timer(); rng = random.Random(5)

            def rand_tails(n):
                for i in range(n):
                    if i % 3 == 0:
                        vals = [rng.randint(1, 1000) for _ in range(3)]
                        yield tuple(rng.choice(vals) for _ in range(7))
                    else:
                        yield tuple(rng.randint(1, 1000) for _ in range(7))
            try:
                n1 = tail_coverage(7, tails_with_sum_at_most(7, 20), leaves)
                n2 = tail_coverage(7, rand_tails(20000), leaves); ok = n1 == 77520
            except (AssertionError, KeyError) as e:
                n1, n2, ok = str(e), 0, False
            res["r7_tails"] = {"exhaustive_sum_bound": 20, "exhaustive_tails": n1, "random_tails": n2}
            checks["tails_r7"] = report("r=7: integer tails (sum <= 20, and random) lie in certified cells with the same verdict",
                                        ok, f"{n1} + {n2} tails ({t()} s)")
    if not a.no_generator and (shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")):
        t = Timer(); exe = build_c("sign_cell_certificates", os.path.join(out, "bin")); ok = True
        for r in range(1, min(a.rmax, 6) + 1):
            got = subprocess.run([exe, str(r), "1", "0", "1"], capture_output=True, check=True).stdout
            with gzip.open(stream_path(r), "rb") as f:
                ok &= got == f.read()
        checks["generator"] = report("the C generator reproduces the streams for r <= 6 byte for byte", ok, f"({t()} s)")
        if a.rmax >= 7:
            t = Timer(); meta = read_json(os.path.join(DATA, "cells", "sign_certificates_r7.json"))
            h = hashlib.sha256(); nbytes = 0
            with subprocess.Popen([exe, "7", "1", "0", "1"], stdout=subprocess.PIPE) as proc:
                for chunk in iter(lambda: proc.stdout.read(1 << 22), b""):
                    h.update(chunk); nbytes += len(chunk)
            good = proc.returncode == 0 and (h.hexdigest(), nbytes) == (meta["uncompressed_sha256"], meta["uncompressed_bytes"])
            res["generator_r7"] = {"uncompressed_sha256": h.hexdigest(), "uncompressed_bytes": nbytes}
            checks["generator_r7"] = report("the C generator reproduces the stream for r = 7 byte for byte (uncompressed)",
                                            good, f"{nbytes} bytes ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "cells_depth7")


if __name__ == "__main__":
    main()
