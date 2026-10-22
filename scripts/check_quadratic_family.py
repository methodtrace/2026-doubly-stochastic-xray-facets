#!/usr/bin/env python3
"""The quadratic-height family G_ab: facet rays exactly for a <= b + 1, of height (b - 1)a + 1.

Claims (thm:height, proved in the paper): for a >= 1, b >= 2 let G_ab be the normalized potential on {0..a + 2b}
with increments (b^a, 0^(b-1), (b-1)a + 1, 1^b). Its tight graph is connected, has total support (a tight permutation
fixes the universal index e = a + b - 1), and G_ab satisfies P_L (witness a) and P_R (witness e); dim Z(G_ab) = 3 if
a <= b + 1 and 4 otherwise. So G_ab is a facet ray iff a <= b + 1, and then a core of height (b - 1)a + 1;
G_{b+1,b} has height b^2 on 3b + 2 points (a = 1 gives the family g_b of check_unbounded_height.py). Consequently, for every k >= 8
some facet ray of hat-P(k) has height at least floor((k - 2)/3)^2 (pad_L^j G_{b+1,b}, b = floor((k - 2)/3),
j = k - 3b - 2).
Checks (exact):
  1. for 2 <= b <= --bmax and 1 <= a <= b + 4, with w = G_ab - iota: v equal to the closed form of the proof, the
     tight graph equal to the stated list (including (0, a + e) iff a <= b + 1), the stated tight permutation fixing
     e, e tight with every index, connected bipartite tight graph, dim Z by rational elimination, minimum increment 0,
     gcd 1, height (b - 1)a + 1, P_L and P_R witnesses, first increment > 0 and last < height (core); for b <= 6 also
     every tight edge allowed (forced matching) and the face dimension 2k - 1 - dim Z by rank (facet iff
     a <= b + 1); for b <= --classify-bmax the certified facet test facet_test.classify agrees (facet iff a <= b + 1);
  2. for 8 <= k <= --kmax: pad_L^j G_{b+1,b} is a ray with total support, connected tight graph and height
     floor((k - 2)/3)^2; facet by facet_test.classify; for k <= 20 also the exact face dimension 2k - 4.
--bmax 12: 121 pairs (88 facet rays, 33 with dim Z = 4).
Default: b <= 8, classify b <= 5, k <= 20 (about 5 s). --full: b <= 12, classify b <= 8, k <= 30 (about 10 s).
"""
import argparse

from common import Timer, finish, out_dir, report

from dsxray import exact_graphs as GR
from dsxray import potential_families as PF
from dsxray.facet_test import classify


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/quadratic_family")
    ap.add_argument("--bmax", type=int, default=8)
    ap.add_argument("--classify-bmax", type=int, default=5)
    ap.add_argument("--kmax", type=int, default=20)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.bmax = max(a.bmax, 12); a.classify_bmax = max(a.classify_bmax, 8); a.kmax = max(a.kmax, 30)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); rows = []; ok = True
    for b in range(2, a.bmax + 1):
        for aa in range(1, b + 5):
            try:
                zd, H, k = PF.quadratic(aa, b, forced_edges=b <= 6)
            except AssertionError:
                ok = False; continue
            row = {"a": aa, "b": b, "k": k, "height": H, "zdim": zd, "facet_ray": zd == 3}
            if b <= a.classify_bmax:
                c = classify(PF.quadratic_potential(aa, b)[1])
                row["certified_facet"] = c["facet"]
                ok &= c["facet"] == (aa <= b + 1) and c["zdim"] == zd and c["total_support"]
            ok &= (zd == 3) == (aa <= b + 1)
            rows.append(row)
    nf = sum(r["facet_ray"] for r in rows)
    res["family"] = {"bmax": a.bmax, "pairs": len(rows), "facet_rays": nf, "rows": rows}
    if a.bmax == 12:
        ok &= (len(rows), nf) == (121, 88)
    checks["family"] = report(f"G_ab, 2 <= b <= {a.bmax}, 1 <= a <= b + 4: v, tight graph, total support, dim Z = 3 iff a <= b + 1",
                              ok, f"{len(rows)} pairs, {nf} facet rays ({t()} s)")

    t = Timer(); ok = True; hb = []
    for k in range(8, a.kmax + 1):
        b = (k - 2) // 3; j = k - (3 * b + 2)
        d, _, _ = PF.quadratic_potential(b + 1, b)
        u = GR.from_increments((0,) * j + tuple(d))
        ray, ts, conn, E = PF.ray_ts_connected(u)
        good = ray and ts and conn and max(d) == b * b and len(u) == k
        good &= classify(u)["facet"]
        if k <= 20:
            good &= GR.dimensions(k, E)["face_dim"] == 2 * k - 4
        ok &= good; hb.append({"k": k, "b": b, "left_paddings": j, "height": b * b})
    res["height_bound"] = hb
    checks["height_bound"] = report(f"facet ray of height floor((k-2)/3)^2 for 8 <= k <= {a.kmax}", ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "quadratic_family")


if __name__ == "__main__":
    main()
