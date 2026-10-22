#!/usr/bin/env python3
"""thm:pad (padding), cor:pad, cor:cores and tab:cores (cores), on the complete facet lists (k <= 12).

  1. The hypotheses of cor:pad on every facet ray (3 <= k <= 12): connected tight graph with total support, and
     the conditions P_L and P_R (the exact check stated after cor:cores and in cor:finer).
  2. thm:pad(a) on every facet ray with k <= --proof-kmax: the min-plus self-convolution and the tight graph of
     padL(u) are those predicted from u.
  3. padL maps every facet ray at k to a facet ray at k+1, k = 5..(KMAX-1) (and k = KMAX with --pad-beyond, which
     classifies the padded rays at k+1 directly).
  4. Cores: every facet ray is padL^a padR^b of a unique core; tab:cores (cores and non-step facets obtained by
     padding, k = 6..12), and the count of cor:cores, |F_k| = sum_j c_j (k - j + 1). The cores are compared
     with data/cores.
  5. Block lifts of the cores (data/cores/block_lifts_kNN.json; deposited data, not cited in the paper): recomputed
     for k <= --block-lifts-kmax and compared with the deposited classification.
Default: k <= 12, proof steps k <= 10, block lifts k <= 11 (about 15 s); --full: proof steps k <= 11, --pad-beyond,
block lifts k <= 12.
"""
import argparse
import os

from common import DATA, Timer, finish, load_facet_rows, load_rays, out_dir, read_json, report, write_json

from dsxray import padding as PD
from dsxray.block_lifts import block_lifts
from dsxray.facet_test import classify

TABLE = {6: (34, 12, 4, 4, 0), 7: (72, 18, 13, 5, 8), 8: (172, 62, 58, 36, 22), 9: (409, 137, 195, 92, 103),
         10: (1014, 368, 622, 290, 332), 11: (2612, 993, 1902, 853, 1049), 12: (7018, 2808, 5748, 2566, 3182)}


def facet_patterns(k):
    if k <= 4:
        return PD.small_facet_patterns(k)
    if k <= 8:
        return {PD.pattern(u) for u in load_rays(k) if classify(list(u))["facet"]}
    return {PD.pattern(r["u"]) for r in load_facet_rows(k) if r["facet"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/padding")
    ap.add_argument("--kmax", type=int, default=12)
    ap.add_argument("--proof-kmax", type=int, default=10)
    ap.add_argument("--pad-beyond", action="store_true")
    ap.add_argument("--write-cores", action="store_true")
    ap.add_argument("--block-lifts-kmax", type=int, default=11)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 12); a.proof_kmax = max(a.proof_kmax, 11); a.pad_beyond = True
        a.block_lifts_kmax = max(a.block_lifts_kmax, 12)
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer(); F = {k: facet_patterns(k) for k in range(3, a.kmax + 1)}
    print(f"   facet patterns loaded ({t()} s)", flush=True)

    t = Timer(); hyp = {}
    for k in range(3, a.kmax + 1):
        bad = 0
        for p in F[k]:
            ts, conn = PD.connected_with_total_support(PD.from_pattern(p))
            bad += not (ts and conn and PD.P_L(p) and PD.P_R(p))
        hyp[k] = bad
    res["hypotheses_failures"] = hyp
    checks["hypotheses"] = report(f"facet rays 3 <= k <= {a.kmax}: connected, total support, P_L and P_R", not any(hyp.values()), f"({t()} s)")

    t = Timer(); bad = 0; n = 0
    for k in range(5, a.proof_kmax + 1):
        for p in F[k]:
            n += 1; bad += not PD.padding_structure_ok(p)
    checks["proof_steps"] = report(f"padding: v' and tight graph as predicted, k <= {a.proof_kmax}", bad == 0, f"{n} facet rays ({t()} s)")

    t = Timer(); pad = {}
    for k in range(5, a.kmax):
        pad[k] = sum(PD.pad_left(p) in F[k + 1] for p in F[k]) == len(F[k])
    if a.pad_beyond:
        k = a.kmax
        pad[k] = all(PD.is_facet_pattern(PD.pad_left(p)) for p in F[k])
    res["padding_maps_facets_to_facets"] = pad
    checks["padding"] = report("padL maps facet rays to facet rays", all(pad.values()), f"k = {min(pad)}..{max(pad)} ({t()} s)")

    t = Timer(); table, cores = PD.core_table(F, 3, a.kmax)
    res["cores"] = table
    ok = all(v["nonunique_core"] == 0 for v in table.values())
    for k in range(6, a.kmax + 1):
        row = tuple(table[k][c] for c in ("ray_facets", "cores", "nonstep_facets", "nonstep_cores", "nonstep_from_padding"))
        ok &= row == TABLE[k]
        csz = {j: len(cores[j]) for j in range(3, k + 1)}
        ok &= len(F[k]) == sum(csz[j] * (k - j + 1) for j in csz)
    dep_ok = True
    for k in range(3, a.kmax + 1):
        p = os.path.join(DATA, "cores", f"cores_k{k:02d}.json")
        if os.path.exists(p):
            dep_ok &= sorted(map(tuple, read_json(p)["cores"])) == cores[k]
        if a.write_cores:
            write_json(os.path.join(out_dir(os.path.join(out, "cores")), f"cores_k{k:02d}.json"),
                       {"k": k, "count": len(cores[k]), "representation": "increment pattern (min 0, gcd 1)",
                        "cores": [list(c) for c in cores[k]]})
    checks["cores"] = report(f"unique cores, core table and |F_k| = sum_j c_j (k-j+1), k <= {a.kmax}", ok, f"({t()} s)")
    checks["cores_match_deposited"] = report("cores equal data/cores", dep_ok)

    t = Timer(); ok = True; bl = {}
    for k in range(5, min(a.block_lifts_kmax, a.kmax) + 1):
        dep = read_json(os.path.join(DATA, "cores", f"block_lifts_k{k:02d}.json"))
        rows = [{"pattern": list(c), "block_lifts": block_lifts(c), "primitive": None, "step": PD.is_step(c)} for c in cores[k]]
        for r in rows:
            r["primitive"] = not r["block_lifts"]
        ok &= rows == dep["cores"]
        bl[k] = {"primitive": sum(r["primitive"] for r in rows), "nonstep_primitive": sum(r["primitive"] and not r["step"] for r in rows)}
        ok &= bl[k]["primitive"] == dep["counts"]["primitive"]
    res["block_lifts"] = bl
    checks["block_lifts"] = report(f"block lifts of the cores recomputed, equal to data/cores/block_lifts, k <= {min(a.block_lifts_kmax, a.kmax)}", ok, f"({t()} s)")
    res["checks"] = checks
    finish(res, out, "padding")


if __name__ == "__main__":
    main()
