#!/usr/bin/env python3
"""Statement (S1) for every real sequence with r = 10 steps: the domain index and the verification of the deposit.

Claim (computer-assisted): (S1) holds for every real sequence E_0..E_10 (see check_s1_certificates.py for the
statement and the method). The certificate tree for r = 10 is split into 441 domains, one per mask prefix of length 6
(radices 1, 1, 3, 3, 7, 7), each an independent certificate tree for its domain; together 123,173,211 terminal
records (C 47,898,083, D 23,959,889, G 26,825,757, Z 24,489,482; 4,664,113,075 bytes uncompressed, 355,486,968 bytes
xz-compressed). Because of their size the domain files are not in this archive: they are in the separate deposit file
s1_depth10_certificates.zip (directory s1_depth10_certificates/, one file <prefix>.txt.xz per domain, stored without
further compression). data/s1/s1_depth10_domains.json is the index: for each domain its prefix, mixed-radix interval,
file name, sizes of the compressed and of the uncompressed stream, and its terminal counts;
data/s1/s1_depth10_certificates.sha256 and data/s1/s1_depth10_trees.sha256 hold the SHA-256 of the 441 compressed
files and of the uncompressed trees (sha256sum format).

Checks:
  1. (always) the index: the 441 prefixes are exactly all mask prefixes of length 6, the intervals are those of the
     prefixes and partition [0, 572834992892439375) (the product of the radices of the 20 levels), and the counts add
     up to the totals;
  2. (--deposit PATH, a zip file or a directory with the 441 files) every file is present; every
     domain tree is verified completely as in check_s1_certificates.py (partition of its domain, every certificate,
     summary counts, uncompressed SHA-256 and terminal counts equal to the index). About 70 core-minutes
     (35 minutes with --workers 2); --domains N verifies only the N smallest domains (a quick test of the deposit).
  3. --make-deposit SRC OUT.zip builds the deposit file from a directory SRC (or several, separated by commas) holding the 441 xz-compressed domain
     trees under any names (each is identified by the SHA-256 of its uncompressed tree; the compressed bytes may
     differ from the index, since they depend on the version of the xz library), writes the checksums of the files
     actually deposited, and checks the deposit as in 2 (--no-verify: presence only).
Default: check 1 (seconds).
"""
import argparse
import hashlib
import json
import lzma
import os
import zipfile
from itertools import product

from common import DATA, Timer, finish, out_dir, pmap, read_json, report

from dsxray import s1_certificates as S1

INDEX = os.path.join(DATA, "s1", "s1_depth10_domains.json")
SHA_FILES = {"compressed_sha256": os.path.join(DATA, "s1", "s1_depth10_certificates.sha256"),
             "uncompressed_sha256": os.path.join(DATA, "s1", "s1_depth10_trees.sha256")}


def load_index():
    """The index with the SHA-256 values of the two checksum files attached to each domain."""
    idx = read_json(INDEX)
    for key, path in SHA_FILES.items():
        with open(path) as f:
            sums = dict(reversed(line.split()) for line in f if line.strip())
        for d in idx["domain_list"]:
            name = d["file"] if key == "compressed_sha256" else d["file"][:-3]
            d[key] = sums[name]
    return idx


def file_sha(f):
    h = hashlib.sha256()
    for chunk in iter(lambda: f.read(1 << 22), b""):
        h.update(chunk)
    return h.hexdigest()


class Deposit:
    def __init__(self, path, folder):
        self.path, self.folder = path, folder
        self.zip = path.endswith(".zip")

    def open(self, name):
        if self.zip:
            return zipfile.ZipFile(self.path).open(f"{self.folder}/{name}")
        return open(os.path.join(self.path, name), "rb")


def _verify_domain(args):
    path, folder, dom = args
    dep = Deposit(path, folder)
    try:
        s = S1.verify_tree(lambda: lzma.open(dep.open(dom["file"])), 10, tuple(dom["prefix"]))
    except AssertionError as e:
        return False, str(e)
    ok = (s["domain_interval"] == dom["interval"] and s["uncompressed_sha256"] == dom["uncompressed_sha256"]
          and s["uncompressed_bytes"] == dom["uncompressed_bytes"] and s["terminal_records"] == dom["terminal_records"]
          and all(s["certificate_types"].get(k, 0) == v for k, v in dom["certificate_types"].items()))
    return ok, s["terminal_records"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/s1_depth10")
    ap.add_argument("--deposit", default="")
    ap.add_argument("--domains", type=int, default=0, help="verify only this many domains (smallest first); 0 = all")
    ap.add_argument("--make-deposit", nargs=2, metavar=("SRC", "OUT_ZIP"))
    ap.add_argument("--no-verify", action="store_true", help="with --deposit or --make-deposit: check presence only")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true", help="no effect without --deposit (the domain files are not in the archive)")
    a = ap.parse_args()
    out = out_dir(a.out); res = {}; checks = {}
    idx = load_index(); doms = idx["domain_list"]

    t = Timer()
    groups, radices, widths = S1.tree_shape(10)
    want = [list(p) for p in product(*[range(1, x + 1) for x in radices[:6]])]
    ok = [d["prefix"] for d in doms] == want and len(doms) == idx["domains"] == 441
    ok &= all(list(S1.interval(tuple(d["prefix"]), radices, widths)) == d["interval"] for d in doms)
    ok &= doms[0]["interval"][0] == 0 and all(x["interval"][1] == y["interval"][0] for x, y in zip(doms, doms[1:]))
    ok &= doms[-1]["interval"][1] == widths[0] == idx["tree_interval"][1] == 572834992892439375
    ok &= sum(d["terminal_records"] for d in doms) == idx["terminal_records"] == 123173211
    ok &= all(sum(d["certificate_types"][k] for d in doms) == idx["certificate_types"][k] for k in "CDGZ")
    ok &= sum(d["uncompressed_bytes"] for d in doms) == idx["uncompressed_bytes"]
    ok &= sum(d["compressed_bytes"] for d in doms) == idx["compressed_bytes"]
    ok &= len({d["file"] for d in doms}) == 441 and len({d["compressed_sha256"] for d in doms}) == 441
    res["index"] = {k: v for k, v in idx.items() if k != "domain_list"}
    checks["index"] = report("r = 10 domain index: 441 domains partition the whole tree, counts add up", ok, f"({t()} s)")

    folder = idx["deposit_directory"]
    if a.make_deposit:
        t = Timer(); src, dst = a.make_deposit
        by_sha = {}
        for p in sorted(os.path.join(s, fn) for s in src.split(",") for fn in os.listdir(s)):
            if os.path.isfile(p) and p.endswith(".txt.xz"):
                with lzma.open(p, "rb") as f:
                    by_sha[file_sha(f)] = p
        missing = [d["file"] for d in doms if d["uncompressed_sha256"] not in by_sha]
        if not missing:
            comp = {}
            for d in doms:
                with open(by_sha[d["uncompressed_sha256"]], "rb") as f:
                    comp[d["file"]] = file_sha(f)
            with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_STORED) as z:
                z.write(INDEX, f"{folder}/s1_depth10_domains.json")
                z.write(SHA_FILES["uncompressed_sha256"], f"{folder}/s1_depth10_trees.sha256")
                z.writestr(f"{folder}/s1_depth10_certificates.sha256",
                           "".join(f"{comp[d['file']]}  {d['file']}\n" for d in doms))
                for d in doms:
                    z.write(by_sha[d["uncompressed_sha256"]], f"{folder}/{d['file']}")
            a.deposit = dst
            res["make_deposit"] = {"files": len(doms),
                                   "compressed_differ_from_index": sum(comp[d["file"]] != d["compressed_sha256"] for d in doms)}
        checks["make_deposit"] = report(f"deposit file written: {dst} (files identified by their uncompressed SHA-256)",
                                        not missing, f"missing {len(missing)} ({t()} s)")
    if a.deposit:
        t = Timer(); dep = Deposit(a.deposit, folder); absent = []; recompressed = []
        for d in doms:
            try:
                with dep.open(d["file"]) as f:
                    if file_sha(f) != d["compressed_sha256"]:
                        recompressed.append(d["file"])
            except (KeyError, FileNotFoundError):
                absent.append(d["file"])
        res["deposit_files"] = {"absent": len(absent), "compressed_differ_from_index": len(recompressed)}
        checks["deposit_files"] = report("every domain file is in the deposit", not absent,
                                         f"{len(absent)} absent; {len(recompressed)} with compressed bytes other than the index "
                                         f"(allowed: identity is checked on the uncompressed tree) ({t()} s)")
        if not absent and not a.no_verify:
            t = Timer()
            sel = sorted(doms, key=lambda d: d["uncompressed_bytes"])[:a.domains] if a.domains else doms
            outs = pmap(_verify_domain, [(a.deposit, folder, d) for d in sel], a.workers, chunksize=1)
            fails = [d["file"] for d, (good, _) in zip(sel, outs) if not good]
            n = sum(x for good, x in outs if good)
            res["deposit"] = {"domains_verified": len(sel), "terminal_records": n, "failures": fails[:20]}
            checks["deposit_domains"] = report(f"{len(sel)} domain trees verified completely (uncompressed SHA-256, "
                                               "partition, every certificate, counts)", not fails,
                                               f"{n} terminal records ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "s1_depth10")


if __name__ == "__main__":
    main()
