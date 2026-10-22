#!/usr/bin/env python3
"""Fingerprint of a paper-tier run, and comparison of two runs or of a run with the deposited data.

  python3 tools/fingerprint.py OUT                    write OUT/fingerprint.json and print a summary
  python3 tools/fingerprint.py --compare A.json B.json   compare two fingerprints (for example two machines)

For every piece folder OUT/pieces/<id>/ the fingerprint holds
  - status: the status recorded in piece.json (timings and memory in piece.json are not part of the fingerprint);
  - results: for each result file <name>.json, the SHA-256 of its canonical form: keys sorted, and the fields that
    describe how a run was made rather than what it computed (paths, worker counts, timings) removed;
  - files: for every other file the piece wrote (rays, facet lists, certificate trees and streams, deposit file), the
    SHA-256 and size of its UNCOMPRESSED content (.gz and .xz files are decompressed; a zip deposit is fingerprinted
    member by member, each member decompressed). Compressed bytes depend on the version of the zlib or xz library
    and are not compared. Logs, compiled programs and resume state are skipped.
Where a regenerated file has a deposited counterpart in data/ (the ray lists and the facet lists), the fingerprint
also records whether its uncompressed content equals the deposited file.

Reproducibility standard: two runs agree when every piece passed in both and every results and files entry is
identical. Floating-point cross-checks (cddlib, the floating-point LP comparison of the cell elimination) report
only pass or fail, so their result files contain no floating-point values.
"""
import argparse
import gzip
import hashlib
import json
import lzma
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VOLATILE = {"args", "out", "workers", "seconds", "wall_seconds", "cpu_seconds", "elapsed", "time", "finished"}
SKIP_DIRS = {"bin", "__pycache__"}
SKIP_FILES = {"piece.json", "log.txt"}


def _canon(o, out):
    if isinstance(o, dict):
        return {k: _canon(v, out) for k, v in o.items() if k not in VOLATILE}
    if isinstance(o, list):
        return [_canon(v, out) for v in o]
    if isinstance(o, str) and out and out in o:
        return o.replace(out, "<OUT>")
    return o


def _stream_sha(f):
    h = hashlib.sha256(); n = 0
    for chunk in iter(lambda: f.read(1 << 22), b""):
        h.update(chunk); n += len(chunk)
    return h.hexdigest(), n


def _open_uncompressed(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    if path.endswith(".xz"):
        return lzma.open(path, "rb")
    return open(path, "rb")


def file_fingerprint(path):
    if path.endswith(".zip"):
        members = {}
        with zipfile.ZipFile(path) as z:
            for name in sorted(z.namelist()):
                if name.endswith("/"):
                    continue
                with z.open(name) as m:
                    f = lzma.open(m) if name.endswith(".xz") else gzip.open(m) if name.endswith(".gz") else m
                    sha, n = _stream_sha(f)
                members[name] = {"sha256_uncompressed": sha, "bytes_uncompressed": n}
        h = hashlib.sha256(json.dumps(members, sort_keys=True).encode()).hexdigest()
        return {"sha256_uncompressed": h, "members": len(members), "member_list": members}
    with _open_uncompressed(path) as f:
        sha, n = _stream_sha(f)
    return {"sha256_uncompressed": sha, "bytes_uncompressed": n}


def deposited_counterpart(rel):
    """The deposited file that a regenerated file must equal (uncompressed), or None."""
    base = os.path.basename(rel)
    stem = base[:-3] if base.endswith(".gz") else base
    if stem.startswith("facets_k") and stem.endswith(".jsonl"):
        p = os.path.join(ROOT, "data", "facets", stem)
        return p if os.path.exists(p) else p + ".gz" if os.path.exists(p + ".gz") else None
    if stem.startswith("rays_k") and stem.endswith(".json"):
        p = os.path.join(ROOT, "data", "rays", stem)
        return p if os.path.exists(p) else None
    return None


def fingerprint(out):
    out = os.path.abspath(out); pieces = {}
    pdir = os.path.join(out, "pieces")
    for pid in sorted(os.listdir(pdir)):
        d = os.path.join(pdir, pid); rec = {"status": None, "results": {}, "files": {}}
        try:
            with open(os.path.join(d, "piece.json")) as f:
                rec["status"] = json.load(f)["status"]
        except (OSError, ValueError, KeyError):
            rec["status"] = "missing"
        for dirpath, dirnames, filenames in os.walk(d):
            dirnames[:] = sorted(x for x in dirnames if x not in SKIP_DIRS)
            for fn in sorted(filenames):
                p = os.path.join(dirpath, fn); rel = os.path.relpath(p, d)
                if rel in SKIP_FILES or "state" in fn or fn.endswith(".part"):
                    continue
                if dirpath == d and fn.endswith(".json"):
                    with open(p) as f:
                        obj = _canon(json.load(f), out)
                    rec["results"][fn] = hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()
                    continue
                fp = file_fingerprint(p)
                dep = deposited_counterpart(rel)
                if dep:
                    fp["equals_deposited"] = file_fingerprint(dep)["sha256_uncompressed"] == fp["sha256_uncompressed"]
                rec["files"][rel] = fp
        pieces[pid] = rec
    return {"run": out, "pieces": pieces}


def compare(a, b):
    diffs = []
    for pid in sorted(set(a["pieces"]) | set(b["pieces"])):
        x, y = a["pieces"].get(pid), b["pieces"].get(pid)
        if x is None or y is None:
            diffs.append(f"{pid}: only in {'B' if x is None else 'A'}"); continue
        if x["status"] != "pass" or y["status"] != "pass":
            diffs.append(f"{pid}: status {x['status']} / {y['status']}")
        for kind in ("results", "files"):
            for name in sorted(set(x[kind]) | set(y[kind])):
                u, v = x[kind].get(name), y[kind].get(name)
                su = u if kind == "results" else (u or {}).get("sha256_uncompressed")
                sv = v if kind == "results" else (v or {}).get("sha256_uncompressed")
                if su != sv:
                    diffs.append(f"{pid}: {kind} {name} differs")
    return diffs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", nargs="?")
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"))
    a = ap.parse_args()
    if a.compare:
        with open(a.compare[0]) as f:
            A = json.load(f)
        with open(a.compare[1]) as f:
            B = json.load(f)
        diffs = compare(A, B)
        n = len(set(A["pieces"]) | set(B["pieces"]))
        for line in diffs:
            print("DIFF", line)
        print(f"{n} pieces compared: {'IDENTICAL' if not diffs else f'{len(diffs)} differences'}")
        return 0 if not diffs else 1
    if not a.out:
        ap.error("give OUT or --compare A B")
    fp = fingerprint(a.out)
    with open(os.path.join(a.out, "fingerprint.json"), "w") as f:
        json.dump(fp, f, indent=1, sort_keys=True)
    P = fp["pieces"]
    failed = [p for p, r in P.items() if r["status"] != "pass"]
    notdep = [f"{p}/{n}" for p, r in P.items() for n, x in r["files"].items() if x.get("equals_deposited") is False]
    print(f"{len(P)} pieces; {len(P) - len(failed)} passed; failed or missing: {', '.join(failed) or 'none'}")
    print(f"regenerated files that differ from the deposited data: {', '.join(notdep) or 'none'}")
    print(f"written: {os.path.join(a.out, 'fingerprint.json')}")
    return 0 if not failed and not notdep else 1


if __name__ == "__main__":
    sys.exit(main())
