#!/usr/bin/env python3
"""Step potentials (staircases): certified sweeps behind the statements about step facets.

For every composition g of k with 1 <= r <= k-2 (staircase u(g), sign +; every signed step potential is of this
form up to an affine term):
  1. certified classification: facet iff total support (conj:T, tab:evidence;
     16,777,168 staircases for 3 <= k <= 24), the three-step rule on every four-part composition, and staircases
     with the same maximizer pattern have the same outcome;
  2. the rigid-corner certificate (conj:E): every staircase with a tight permutation but without total support
     has a level with a unique attaining pair whose corner lies in no tight permutation (4,582,776 for k <= 24); and
     for the first block start with a single tight neighbor E_q, 0 < q < r, a unique-attainer partner exists;
  3. the criterion of thm:S4 (S4 in the code) never admits a non-facet (764 facet compositions through g or g*
     for 3 <= k <= 12, 308 of them outside the families of thm:S), g -> g* is an involution that reflects the
     tight graph, the tight neighbors of block starts are the attaining ends (lem:starts), and the tight pairs on the
     anti-diagonal M(s) are the corners of the attaining block pairs (lem:last).
Staircases with a nonincreasing tail (thm:dec) are checked by check_nonincreasing_tails.py.
Defaults: sweep and rigid corner k <= 13, part 3 k <= 12. --full: sweep and rigid corner k <= 24 (100 to 150
core-hours, dominated by k = 23, 24).

Running in pieces: --k K --r R1:R2 runs parts 1 and 2 only for the staircases of k with R1 <= r <= R2 steps
(--no-small skips part 3), and --merge DIR [DIR ...] combines the result files of such pieces (each DIR holds a
staircases.json) into the full result and applies every check. The pieces must cover every (k, r) with 3 <= k <= KMAX
and 1 <= r <= k - 2 exactly once; the maximizer patterns are keyed by r, so the counts of the pieces add up exactly.
A run without --k gives the same result as the pieces merged.
"""
import argparse
from collections import Counter

import glob
import json
import os

from common import DATA, Timer, finish, out_dir, pimap, read_json, report

from dsxray import level_graph as LG
from dsxray import staircases as ST
from dsxray.facet_test import classify

COUNTS = {k: 2 ** (k - 1) - 2 for k in range(3, 25)}
RIGID_TOTAL_K24 = 4582776


def _task(g):
    r = ST.classify_staircase(g)
    return {key: r[key] for key in ("g", "r", "facet", "total_support", "ray", "cert", "triple_rule") if key in r}


def _rigid(g):
    return LG.rigid_corner(g), LG.first_single_start_partner(g)


def fam_thmS(J, k):
    J = sorted(J); r = len(J)
    return r <= 2 or J == list(range(J[0], J[0] + r))


SWEEP_SUM = ("items", "facets", "total_support", "facet_equals_total_support", "rays", "three_step",
             "three_step_rule_agrees", "maximizer_patterns", "patterns_with_both_outcomes")


def sweep_part(k, rs, workers):
    """Parts 1 and 2 for the staircases of k with r in rs: (sweep summary, rigid-corner counts)."""
    gen = (g for r in rs for g in ST.staircases_r(k, r))
    s = ST.summarize_classifications(pimap(_task, gen, workers, chunksize=256))
    s = {key: s.get(key, 0) for key in SWEEP_SUM} | {"nonfacet_certificates": s["nonfacet_certificates"]}
    C = Counter(); SS = Counter()
    for c, f in pimap(_rigid, (g for r in rs for g in ST.staircases_r(k, r)), workers, chunksize=256):
        C[c] += 1; SS[str(f)] += 1
    return s, {"counts": dict(sorted(C.items())), "first_single_start": dict(sorted(SS.items()))}


def add_sweep(a, b):
    out = {key: a.get(key, 0) + b.get(key, 0) for key in SWEEP_SUM}
    nf = Counter(a.get("nonfacet_certificates", {})); nf.update(b.get("nonfacet_certificates", {}))
    out["nonfacet_certificates"] = dict(sorted(nf.items()))
    return out


def add_rigid(a, b):
    out = {}
    for key in ("counts", "first_single_start"):
        c = Counter(a.get(key, {})); c.update(b.get(key, {})); out[key] = dict(sorted(c.items()))
    return out


def sweep_checks(res, kmax, checks):
    tot = 0; certified = 0
    for k in range(3, kmax + 1):
        s = res["sweep"][k]
        ok = s["items"] == COUNTS[k] and s["facet_equals_total_support"] == s["items"] and s["patterns_with_both_outcomes"] == 0
        ok = ok and s.get("three_step", 0) == s.get("three_step_rule_agrees", 0)
        ok = ok and set(s["nonfacet_certificates"]) <= {"empty", "coords"}
        tot += s["items"]
        checks[f"sweep_k{k}"] = report(f"k={k} facet iff total support on all staircases",
                                       ok, f"{s['items']} staircases, {s['facets']} facets")
        C = res["rigid_corner"][k]["counts"]; SS = res["rigid_corner"][k]["first_single_start"]
        certified += C.get("certified", 0)
        checks[f"rigid_corner_k{k}"] = report(f"k={k} rigid-corner certificate", C.get("FAIL", 0) == 0 and SS.get("False", 0) == 0,
                                              str(C))
    res["sweep_total"] = tot
    res["rigid_corner_certified_total"] = certified
    if kmax >= 24:
        checks["rigid_corner_total"] = report("rigid-corner certified total k <= 24", certified == RIGID_TOTAL_K24, str(certified))


def part3(small_kmax, res, checks):
    s4 = s4_outside = facets = 0; bad = 0
    for k in range(3, small_kmax + 1):
        for g in ST.staircases(k):
            f = classify(ST.u_of(g))["facet"]; facets += f
            gs = ST.star(g)
            bad += (ST.star(gs) != g) or (classify(ST.u_of(gs))["facet"] != f) or not LG.gstar_tight_graph_reflects(g)
            bad += (not LG.start_neighbors_are_attaining_ends(g)) or (not LG.last_antidiagonal_corners(g))
            if ST.S4(g) or ST.S4(gs):
                s4 += 1; bad += (not f)
                J = sorted(ST.J_of(g)); Jc = [j for j in range(1, k) if j not in J]
                s4_outside += (not fam_thmS(J, k) and not fam_thmS(Jc, k))
    res["S4"] = {"k_range": [3, small_kmax], "facet_compositions": facets, "certified_by_S4": s4,
                 "certified_outside_families": s4_outside}
    checks["S4_and_lemmas"] = report(f"S4 sufficient, g -> g* involution and block lemmas, k <= {small_kmax}", bad == 0,
                                     f"S4 certifies {s4} of {facets} facet compositions, {s4_outside} outside the families")
    if small_kmax == 12:
        checks["S4_counts"] = report("S4 counts for 3 <= k <= 12", (s4, s4_outside, facets) == (764, 308, 2811))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/staircases")
    ap.add_argument("--kmax", type=int, default=13)
    ap.add_argument("--small-kmax", type=int, default=12)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--k", type=int, default=0, help="run parts 1 and 2 for this k only (a piece)")
    ap.add_argument("--r", default="", help="R1:R2, the step counts of the piece (default: all)")
    ap.add_argument("--no-small", action="store_true", help="skip part 3")
    ap.add_argument("--merge", nargs="+", default=None, metavar="DIR", help="combine the results of pieces")
    a = ap.parse_args()
    if a.full:
        a.kmax = max(a.kmax, 24)
    res = {"sweep": {}, "rigid_corner": {}}; checks = {}

    if a.merge is not None:
        parts = []
        for d in a.merge:
            for p in sorted(glob.glob(os.path.join(d, "staircases.json"))):
                with open(p) as f:
                    parts.append(json.load(f))
        seen = Counter(); kmax = 0; s4 = None
        for P in parts:
            if "S4" in P:
                s4 = P["S4"]
            for piece in P.get("pieces", []):
                k = piece["k"]; kmax = max(kmax, k)
                for r in range(piece["r"][0], piece["r"][1] + 1):
                    seen[(k, r)] += 1
                res["sweep"][k] = add_sweep(res["sweep"].get(k, {}), piece["sweep"])
                res["rigid_corner"][k] = add_rigid(res["rigid_corner"].get(k, {}), piece["rigid_corner"])
        want = {(k, r) for k in range(3, kmax + 1) for r in range(1, k - 1)}
        cover_ok = set(seen) == want and all(v == 1 for v in seen.values())
        checks["pieces_cover"] = report(f"the pieces cover every (k, r) with 3 <= k <= {kmax} exactly once", cover_ok,
                                        f"{len(seen)} of {len(want)}; repeated {sum(1 for v in seen.values() if v > 1)}")
        res["sweep"] = dict(sorted(res["sweep"].items())); res["rigid_corner"] = dict(sorted(res["rigid_corner"].items()))
        if cover_ok:
            sweep_checks(res, kmax, checks)
            dep = read_json(os.path.join(DATA, "summary", "staircases.json"))
            diff = []
            for k in range(3, kmax + 1):
                want_s = dep["staircases"]["per_k"].get(str(k)); want_r = dep["rigid_corner"]["per_k"].get(str(k))
                if want_s is not None and {key: want_s.get(key, 0) for key in res["sweep"][k]} != res["sweep"][k]:
                    diff.append(f"sweep k={k}")
                if want_r is not None and want_r != res["rigid_corner"][k]:
                    diff.append(f"rigid k={k}")
            checks["deposited_summary"] = report("every per-k count equals data/summary/staircases.json", not diff,
                                                 ", ".join(diff[:10]))
        if s4 is not None:
            res["S4"] = s4
            checks["S4_counts"] = report("S4 counts for 3 <= k <= 12", (s4["certified_by_S4"], s4["certified_outside_families"],
                                                                         s4["facet_compositions"]) == (764, 308, 2811))
        res["checks"] = checks
        finish(res, out_dir(a.out), "staircases")

    if a.k:
        k = a.k
        r1, r2 = (int(x) for x in a.r.split(":")) if a.r else (1, k - 2)
        t = Timer(); s, rc = sweep_part(k, range(r1, r2 + 1), a.workers)
        res = {"pieces": [{"k": k, "r": [r1, r2], "sweep": s, "rigid_corner": rc}]}
        checks["piece"] = report(f"k={k}, r={r1}..{r2}: facet iff total support, rigid corner",
                                 s["facet_equals_total_support"] == s["items"] and s["patterns_with_both_outcomes"] == 0
                                 and set(s["nonfacet_certificates"]) <= {"empty", "coords"}
                                 and rc["counts"].get("FAIL", 0) == 0 and rc["first_single_start"].get("False", 0) == 0,
                                 f"{s['items']} staircases ({t()} s)")
        if not a.no_small:
            part3(a.small_kmax, res, checks)
        res["checks"] = checks
        finish(res, out_dir(a.out), "staircases")

    pieces = []
    for k in range(3, a.kmax + 1):
        s, rc = sweep_part(k, range(1, k - 1), a.workers)
        res["sweep"][k] = s; res["rigid_corner"][k] = rc
        pieces.append({"k": k, "r": [1, k - 2], "sweep": s, "rigid_corner": rc})
    sweep_checks(res, a.kmax, checks)
    if not a.no_small:
        part3(a.small_kmax, res, checks)
    res["pieces"] = pieces
    res["checks"] = checks
    finish(res, out_dir(a.out), "staircases")


if __name__ == "__main__":
    main()
