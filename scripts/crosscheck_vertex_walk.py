#!/usr/bin/env python3
"""Independent recomputation of the ray lists by walking the vertices of Q = sum_t Q_t (shares no code with the
facet walk), and floating-point vertex-facet enumerations for small k.

  1. Exact version (integer double description, no floating point) for 5 <= k <= --kmax: number of vertices of Q and
     the ray set, compared with data/summary/vertex_walk.json and data/rays.
  2. cddlib version (optional dependency pycddlib; rays verified exactly) for 5 <= k <= --cdd-kmax.
  3. cddlib facets of hat-P(k) from all profiles, every row verified exactly, for k <= --dd-kmax (at most 8),
     compared with the exact facet lists (counts 13, 23, 45, 85, 187 for k = 4..8).
Defaults: exact k <= 9, cddlib k <= 9 and profiles k <= 7 when pycddlib is installed. --full: exact k <= 11
(about 45 minutes at k = 11), cddlib k <= 10, profiles k <= 8. The exact walk at k = 12 (--kmax 12) is estimated
at 5 to 10 core-hours and the cddlib walk at k = 12 (--cdd-kmax 12) took about 3 hours on two workers.
--kmin and --cdd-kmin start the exact and the cddlib walk at a given k (to run one k as a piece); --require-cdd
fails instead of skipping the cddlib parts when pycddlib is missing.
"""
import argparse
import os

from common import DATA, Timer, finish, load_rays, out_dir, read_json, report

from dsxray import small_k_fan as SK
from dsxray import vertex_walk as VW

FACETS = {4: 13, 5: 23, 6: 45, 7: 85, 8: 187}


def have_cdd():
    try:
        import cdd  # noqa: F401
        return True
    except ImportError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/vertex_walk")
    ap.add_argument("--kmax", type=int, default=9)
    ap.add_argument("--cdd-kmax", type=int, default=9)
    ap.add_argument("--dd-kmax", type=int, default=7)
    ap.add_argument("--kmin", type=int, default=5, help="first k of the exact walk (to run one k as a piece)")
    ap.add_argument("--cdd-kmin", type=int, default=5, help="first k of the cddlib walk (to run one k as a piece)")
    ap.add_argument("--require-cdd", action="store_true", help="fail instead of skipping when pycddlib is missing")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 11); a.cdd_kmax = max(a.cdd_kmax, 10); a.dd_kmax = max(a.dd_kmax, 8)
    table = {int(k): v for k, v in read_json(os.path.join(DATA, "summary", "vertex_walk.json"))["table"].items()}
    res = {"exact": {}, "cdd": {}, "double_description": {}}; checks = {}
    for k in range(max(5, a.kmin), a.kmax + 1):
        t = Timer(); nv, rays = VW.walk(k, workers=a.workers if k >= 10 else 1)
        same = rays == set(load_rays(k))
        res["exact"][k] = {"vertices": nv, "rays": len(rays), "identical_to_facet_walk": same}
        checks[f"exact_k{k}"] = report(f"k={k} exact vertex walk", same and nv == table[k]["vertices_of_Q"],
                                       f"{nv} vertices, {len(rays)} rays ({t()} s)")
    cdd = have_cdd()
    if not cdd:
        print("pycddlib is not installed: skipping the cddlib checks", flush=True)
        if a.require_cdd:
            checks["pycddlib"] = report("pycddlib is installed (required by --require-cdd)", False)
    else:
        from dsxray import vertex_walk_cdd as VC
        for k in range(max(5, a.cdd_kmin), a.cdd_kmax + 1):
            t = Timer(); nv, rays = VC.walk(k, workers=a.workers if k >= 10 else 1)
            same = rays == set(load_rays(k))
            res["cdd"][k] = {"vertices": nv, "rays": len(rays), "identical_to_facet_walk": same}
            checks[f"cdd_k{k}"] = report(f"k={k} cddlib vertex walk", same and nv == table[k]["vertices_of_Q"],
                                         f"{nv} vertices, {len(rays)} rays ({t()} s)")
        for k in range(4, min(a.dd_kmax, 8) + 1):
            t = Timer(); H = SK.cdd_facets_verified(k); E, _ = SK.exact_facets(k)
            res["double_description"][k] = len(H)
            checks[f"dd_k{k}"] = report(f"k={k} cddlib facets, rows verified exactly", H == E and len(H) == FACETS[k],
                                        f"{len(H)} facets ({t()} s)")
    res["pycddlib_available"] = cdd
    res["checks"] = checks
    finish(res, out_dir(a.out), "vertex_walk")


if __name__ == "__main__":
    main()
