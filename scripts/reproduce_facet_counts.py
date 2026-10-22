#!/usr/bin/env python3
"""Facet counts of hat-P(k), k = 2..12, from the deposited ray lists, with every claim rechecked.

  1. Classify every deposited ray (data/rays) with the certified facet test; add the 2k-1 coordinate inequalities;
     check the facet counts 2, 6, 13, 23, 45, 85, 187, 426, 1033, 2633, 7041, that every listed ray spans a ray of the
     fan (dim Z = 3), that a ray gives a facet exactly when its tight graph has total support (all 113,842 rays for
     6 <= k <= 12), that distinct facet rays give distinct facets and none equals a coordinate facet, and that the
     outcomes and certificates equal the deposited facet files (data/facets).
  2. For k <= --walk-kmax recompute the ray list by the facet walk (exact, provably complete) and compare it with the
     deposited list; the walk must find every ridge from both sides.
  3. For k <= --profiles-kmax (at most 8) compute the exact facet list from all permutation profiles (line
     enumeration of the fan plus an exact audit of every candidate inequality) and compare.
  4. For --line-kmax >= 9 recompute the deposited ray list by line enumeration (17,500 and 342,087 candidate lines
     for k = 9, 10, as in sec:s_census-cross).

Defaults (quick): classification k <= 12, walk k <= 9, profiles k <= 8, no line enumeration beyond k = 8.
--full: classification k <= 12, walk k <= 11, line enumeration k <= 9. The walk at k = 12 (--walk-kmax 12) takes
about 9 core-hours and the line enumeration at k = 10 (--line-kmax 10) about 8 core-hours. --walk-kmin and
--line-kmin start the walk and the line enumeration at a given k, so that one k can run as a separate piece.
"""
import argparse
import os
import sys
from multiprocessing import Pool

from common import Timer, facets_path, finish, load_facet_rows, load_rays, out_dir, report, write_jsonl, write_rays

from dsxray import fan_walk as FW
from dsxray import small_k_fan as SK
from dsxray.facet_test import classify, classify_coordinate

FACETS = {2: 2, 3: 6, 4: 13, 5: 23, 6: 45, 7: 85, 8: 187, 9: 426, 10: 1033, 11: 2633, 12: 7041}
RAYS = {5: 16, 6: 42, 7: 110, 8: 338, 9: 1146, 10: 4370, 11: 18515, 12: 89321}
CANDIDATE_LINES = {9: 17500, 10: 342087}
FIELDS = ("u", "ray", "zdim", "facet", "cert", "total_support", "lost_edges", "empty_antidiagonals", "face_dim")


def _classify(u):
    r = classify(list(u)); r["u"] = list(u)
    return {f: r[f] for f in FIELDS if f in r}


def classify_all(rays, workers):
    if workers > 1:
        with Pool(workers) as pool:
            return pool.map(_classify, rays, chunksize=64)
    return [_classify(u) for u in rays]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/facet_counts")
    ap.add_argument("--kmax", type=int, default=12, help="classify deposited rays for 5 <= k <= KMAX (<= 12)")
    ap.add_argument("--walk-kmax", type=int, default=9, help="recompute the facet walk for 5 <= k <= this")
    ap.add_argument("--profiles-kmax", type=int, default=8, help="exact facet lists from profiles, 2 <= k <= this (<= 8)")
    ap.add_argument("--walk-kmin", type=int, default=5, help="first k of the facet walk (to run one k as a piece)")
    ap.add_argument("--line-kmin", type=int, default=9, help="first k of the line enumeration (to run one k as a piece)")
    ap.add_argument("--line-kmax", type=int, default=8, help="line enumeration for 9 <= k <= this")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--write-facets", action="store_true", help="write the classified rows to OUT/facets")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 12); a.walk_kmax = max(a.walk_kmax, 11); a.line_kmax = max(a.line_kmax, 9)
    out = out_dir(a.out); checks = {}; res = {"per_k": {}}

    # 1. classification of deposited rays
    for k in range(5, a.kmax + 1):
        t = Timer(); rays = load_rays(k)
        rows = classify_all(rays, a.workers)
        coords = [classify_coordinate(k, tt) for tt in range(2 * k - 1)]
        fac = [r for r in rows if r["facet"]]
        nc = sum(c["facet"] for c in coords)
        certs = {}
        for r in rows:
            certs[r["cert"]] = certs.get(r["cert"], 0) + 1
        ineq = {SK.ray_inequality(r["u"]) for r in fac}
        coord_ineq = {SK.coordinate_inequality(k, tt) for tt in range(2 * k - 1)}
        rec = {"rays": len(rays), "ray_facets": len(fac), "coordinate_facets": nc, "facets": len(fac) + nc,
               "certificates": dict(sorted(certs.items())),
               "facet_equals_total_support": sum(r["facet"] == bool(r.get("total_support")) for r in rows),
               "all_rays_have_dim_z_3": all(r["ray"] for r in rows),
               "distinct_facet_inequalities": len(ineq), "facet_rays_equal_to_a_coordinate_facet": len(ineq & coord_ineq)}
        res["per_k"][k] = rec
        checks[f"k{k}_ray_count"] = report(f"k={k} ray count", len(rays) == RAYS[k], f"{len(rays)}")
        checks[f"k{k}_facet_count"] = report(f"k={k} facet count", rec["facets"] == FACETS[k],
                                              f"{len(fac)} ray facets + {nc} coordinate = {rec['facets']}")
        checks[f"k{k}_rays_are_rays"] = report(f"k={k} every listed ray has dim Z = 3", rec["all_rays_have_dim_z_3"])
        checks[f"k{k}_facet_iff_total_support"] = report(
            f"k={k} facet iff total support", rec["facet_equals_total_support"] == len(rays),
            f"{rec['facet_equals_total_support']}/{len(rays)}")
        checks[f"k{k}_distinct"] = report(f"k={k} distinct facet inequalities, none a coordinate facet",
                                          len(ineq) == len(fac) and not (ineq & coord_ineq))
        dep = facets_path(k)
        if os.path.exists(dep):
            old = {tuple(r["u"]): r for r in load_facet_rows(k)}
            same = len(old) == len(rows) and all(old.get(tuple(r["u"])) == r for r in rows)
            checks[f"k{k}_matches_deposited_facets"] = report(f"k={k} outcomes and certificates equal the deposited file", same)
        if a.write_facets:
            fd = out_dir(os.path.join(out, "facets"))
            write_jsonl(os.path.join(fd, f"facets_k{k:02d}.jsonl" + (".gz" if k >= 12 else "")), rows, compress=k >= 12)
        rec["coordinates"] = [{"t": c["t"], "facet": c["facet"], "cert": c["cert"]} for c in coords]
        print(f"   k={k}: {t()} s", flush=True)

    if a.kmax >= 12:
        ratios = [RAYS[k + 1] / RAYS[k] for k in range(5, 12)]
        res["ray_count_ratios"] = [float(f"{x:.3f}") for x in ratios]
        checks["ray_count_ratios"] = report("ratio of consecutive ray counts goes from about 2.6 (k = 5, 6, 7) to 4.8, "
                                            "increasing from k = 6 on",
                                            f"{ratios[0]:.1f}" == f"{ratios[1]:.1f}" == "2.6" and f"{ratios[-1]:.1f}" == "4.8"
                                            and all(x < y for x, y in zip(ratios[1:], ratios[2:])), str(res["ray_count_ratios"]))
        res["rays_6_to_12"] = sum(RAYS[k] for k in range(6, 13))
        checks["rays_6_to_12"] = report("113,842 rays for 6 <= k <= 12", res["rays_6_to_12"] == 113842)

    # 2. facet walk
    for k in range(max(5, a.walk_kmin), a.walk_kmax + 1):
        t = Timer()
        if k <= 11:
            rays, st = FW.fan_rays_parallel(k, a.workers) if a.workers > 1 and k >= 10 else FW.fan_rays(k)
            sym = st.get("adjacency_symmetric", False)
        else:
            state = os.path.join(out_dir(os.path.join(out, "walk")), f"walk_state_k{k:02d}.jsonl")
            rays, st = FW.fan_rays_resumable(k, a.workers, state, walk_rank=60, log=sys.stderr)
            sym = st["one_sided_ridges"] == 0
        same = set(rays) == set(load_rays(k))
        res.setdefault("walk", {})[k] = {"rays": len(rays), "orbit_representatives": st.get("orbit_reps"),
                                         "ridge_incidences_on_representatives": st.get("ridge_incidences_reps"),
                                         "every_ridge_found_from_both_sides": sym}
        checks[f"k{k}_walk"] = report(f"k={k} facet walk reproduces the deposited rays", same and sym,
                                      f"{len(rays)} rays, {t()} s")
        if k >= 10:
            write_rays(os.path.join(out_dir(os.path.join(out, "walk")), f"rays_k{k:02d}.json"), k, rays)

    # 3. exact facet lists from profiles
    for k in range(2, min(a.profiles_kmax, 8) + 1):
        t = Timer(); H, rays = SK.exact_facets(k)
        ok = len(H) == FACETS[k]
        if k >= 5:
            ok = ok and {SK.normalize_ray(u) for u in rays} == set(load_rays(k))
            if k <= a.kmax:
                rows = classify_all(load_rays(k), 1)
                mine = {SK.ray_inequality(r["u"]) for r in rows if r["facet"]}
                mine |= {SK.coordinate_inequality(k, tt) for tt in range(2 * k - 1) if classify_coordinate(k, tt)["facet"]}
                ok = ok and mine == H
        res.setdefault("profiles", {})[k] = {"facets": len(H), "rays_by_line_enumeration": len(rays)}
        checks[f"k{k}_profiles"] = report(f"k={k} exact facet list from all profiles", ok, f"{len(H)} facets, {t()} s")

    # 4. line enumeration beyond k = 8
    for k in range(max(9, a.line_kmin), a.line_kmax + 1):
        t = Timer()
        d, groups, C = SK.fan_setup(k)
        tasks = SK.prefixes(len(C), d)
        with Pool(a.workers) as pool:
            parts = pool.starmap(SK.cross_directions, [(C, d, p) for p in tasks], chunksize=4)
        dirs = set().union(*parts)
        rays = SK.rays_from_directions(k, dirs, groups, d)
        same = {SK.normalize_ray(u) for u in rays} == set(load_rays(k))
        res.setdefault("line_enumeration", {})[k] = {"candidate_lines": len(dirs), "rays": len(rays)}
        same = same and len(dirs) == CANDIDATE_LINES.get(k, len(dirs))
        checks[f"k{k}_line_enumeration"] = report(f"k={k} line enumeration reproduces the deposited rays "
                                                  f"(candidate lines as in sec:s_census-cross)", same,
                                                  f"{len(dirs)} candidate lines, {len(rays)} rays, {t()} s")
    res["checks"] = checks
    finish(res, out, "facet_counts")


if __name__ == "__main__":
    main()
