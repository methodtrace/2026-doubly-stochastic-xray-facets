"""Shared helpers for the scripts: paths, deposited data, output files."""
import gzip
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "src"))


def out_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def _jsonable(o):
    if hasattr(o, "item"):
        return o.item()
    if isinstance(o, (set, frozenset, tuple)):
        return sorted(o) if isinstance(o, (set, frozenset)) else list(o)
    return str(o)


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=_jsonable)
        f.write("\n")


def read_json(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return json.load(f)


def write_rays(path, k, rays):
    """Ray list, one ray per line."""
    rays = sorted(tuple(u) for u in rays)
    with open(path, "w") as f:
        f.write('{"k": %d, "count": %d, "normalization": "u(0) = u(k-1) = 0, primitive integer vector",\n "rays": [\n' % (k, len(rays)))
        f.write(",\n".join("  " + json.dumps(list(u)) for u in rays))
        f.write("\n ]}\n")


def rays_path(k):
    return os.path.join(DATA, "rays", f"rays_k{k:02d}.json")


def load_rays(k):
    return [tuple(u) for u in read_json(rays_path(k))["rays"]]


def facets_path(k):
    p = os.path.join(DATA, "facets", f"facets_k{k:02d}.jsonl")
    return p if os.path.exists(p) else p + ".gz"


def load_facet_rows(k, path=None):
    path = path or facets_path(k)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]


def write_jsonl(path, rows, compress=False):
    opener = (lambda p: gzip.GzipFile(p, "wb", mtime=0)) if compress else (lambda p: open(p, "wb"))
    with opener(path) as f:
        for r in rows:
            f.write((json.dumps(r) + "\n").encode())


class Timer:
    def __init__(self):
        self.t = time.time()

    def __call__(self):
        return float(f"{time.time() - self.t:.1f}")


def report(name, ok, detail=""):
    ok = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    return ok


def finish(results, out, name):
    """Write results to OUT/name.json and exit 0 iff every check passed."""
    ok = all(results.get("checks", {}).values()) if results.get("checks") else True
    results["all_checks_passed"] = ok
    write_json(os.path.join(out_dir(out), name + ".json"), results)
    print(("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED") + f" -> {os.path.join(out, name + '.json')}", flush=True)
    sys.exit(0 if ok else 1)


def pmap(fn, items, workers=1, chunksize=64):
    """Ordered parallel map (multiprocessing) or plain map for one worker."""
    items = list(items)
    if workers > 1 and len(items) > chunksize:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            return pool.map(fn, items, chunksize=chunksize)
    return [fn(x) for x in items]


def pimap(fn, items, workers=1, chunksize=256):
    """Ordered parallel map that yields results as they arrive (results are not collected in memory)."""
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            yield from pool.imap(fn, items, chunksize=chunksize)
    else:
        yield from map(fn, items)


def build_c(name, workdir):
    """Compile tools/c/<name>.c into workdir (cc -O2) and return the path of the executable."""
    import shutil
    import subprocess
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if not cc:
        raise SystemExit(f"a C compiler (cc) is needed for tools/c/{name}.c")
    exe = os.path.join(out_dir(workdir), name)
    subprocess.run([cc, "-O2", "-o", exe, os.path.join(ROOT, "tools", "c", name + ".c")], check=True)
    return exe


def read_jsonl(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]
