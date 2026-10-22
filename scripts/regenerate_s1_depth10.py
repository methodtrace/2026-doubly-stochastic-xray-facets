#!/usr/bin/env python3
"""Statement (S1) for r = 10: regenerate the 441 domain trees from the C generator and verify each one completely.

For each selected domain of data/s1/s1_depth10_domains.json the generator tools/c/s1_certificate_tree.c writes the
certificate tree of the domain; the uncompressed output must have the SHA-256 and size recorded in
data/s1/s1_depth10_trees.sha256 and in the index. The tree is written xz-compressed (preset 3, the setting of the
deposited files) to OUT/trees/<prefix>.txt.xz and then verified completely by the independent Python checker
(dsxray.s1_certificates.verify_tree: partition of the domain interval, every certificate, summary counts, terminal
counts equal to the index). The SHA-256 of the compressed file is reported but is not a check: compressed bytes depend
on the version of the xz library, uncompressed bytes do not.

The domains are taken in the order of the index; --domains A:B selects domains A..B-1 of that order, so the domains
can be regenerated in pieces (for example 0:100, 100:200, ...). The trees written by all pieces form the deposit file
(check_s1_depth10.py --make-deposit OUT/trees deposit.zip).
Default: the 20 smallest domains (seconds). All 441 domains: a few core-hours; peak memory is that of the checker.
"""
import argparse
import hashlib
import lzma
import os
import subprocess

from common import DATA, Timer, build_c, finish, out_dir, pmap, read_json, report

from dsxray import s1_certificates as S1


def load_index():
    idx = read_json(os.path.join(DATA, "s1", "s1_depth10_domains.json"))
    with open(os.path.join(DATA, "s1", "s1_depth10_trees.sha256")) as f:
        sums = dict(reversed(line.split()) for line in f if line.strip())
    for d in idx["domain_list"]:
        d["uncompressed_sha256"] = sums[d["file"][:-3]]
    return idx


def _regenerate(args):
    exe, dom, folder = args
    path = os.path.join(folder, dom["file"])
    prefix = ",".join(str(x) for x in dom["prefix"])
    h = hashlib.sha256(); n = 0
    with lzma.open(path + ".part", "wb", preset=3) as z, \
            subprocess.Popen([exe, "S", "10", "-1", "0", "1", "1", "1", "0", prefix], stdout=subprocess.PIPE) as proc:
        for chunk in iter(lambda: proc.stdout.read(1 << 22), b""):
            h.update(chunk); n += len(chunk); z.write(chunk)
    ok = proc.returncode == 0 and h.hexdigest() == dom["uncompressed_sha256"] and n == dom["uncompressed_bytes"]
    os.replace(path + ".part", path)
    try:
        s = S1.verify_tree(lambda: lzma.open(path, "rb"), 10, tuple(dom["prefix"]))
        ok = ok and (s["domain_interval"] == dom["interval"] and s["uncompressed_sha256"] == dom["uncompressed_sha256"]
                     and s["terminal_records"] == dom["terminal_records"]
                     and all(s["certificate_types"].get(k, 0) == v for k, v in dom["certificate_types"].items()))
    except AssertionError as e:
        ok = False; s = {"error": str(e)}
    with open(path, "rb") as f:
        csha = hashlib.file_digest(f, "sha256").hexdigest()
    return {"file": dom["file"], "ok": ok, "uncompressed_sha256": h.hexdigest(), "uncompressed_bytes": n,
            "terminal_records": s.get("terminal_records"), "compressed_sha256": csha}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/s1_regenerate")
    ap.add_argument("--domains", default="", help="A:B, domains A..B-1 in the order of the index (default: the 20 smallest)")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    out = out_dir(a.out); folder = out_dir(os.path.join(out, "trees")); checks = {}
    idx = load_index(); doms = idx["domain_list"]
    if a.domains:
        lo, hi = (int(x) for x in a.domains.split(":"))
        sel = doms[lo:hi]
    else:
        lo, hi = None, None
        sel = sorted(doms, key=lambda d: d["uncompressed_bytes"])[:20]
    exe = build_c("s1_certificate_tree", os.path.join(out, "bin"))
    t = Timer()
    outs = pmap(_regenerate, [(exe, d, folder) for d in sorted(sel, key=lambda d: -d["uncompressed_bytes"])],
                a.workers, chunksize=1)
    outs = sorted(outs, key=lambda o: o["file"])
    fails = [o["file"] for o in outs if not o["ok"]]
    want_c = {}
    with open(os.path.join(DATA, "s1", "s1_depth10_certificates.sha256")) as f:
        for line in f:
            if line.strip():
                hsh, name = line.split(); want_c[name] = hsh
    recompressed = [o["file"] for o in outs if o["compressed_sha256"] and o["compressed_sha256"] != want_c[o["file"]]]
    res = {"domains": [lo, hi] if a.domains else "20 smallest", "selected": len(sel),
           "terminal_records": sum(o["terminal_records"] or 0 for o in outs),
           "uncompressed_bytes": sum(o["uncompressed_bytes"] for o in outs),
           "compressed_bytes_differ_from_index": len(recompressed), "trees": outs}
    checks["regenerate"] = report(f"{len(sel)} domain trees regenerated byte for byte (uncompressed) and verified completely",
                                  not fails, f"{res['terminal_records']} terminal records, failures {fails[:5]} ({t()} s)")
    finish({"results": res, "checks": checks}, out, "s1_regenerate")


if __name__ == "__main__":
    main()
