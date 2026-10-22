"""Statement E-flat for all tails with r steps: exact verification of a complete sign-cell certificate stream.

This is a second route to thm:Eflat6 (app:cells), independent of the Fourier-Motzkin enumeration of
cell_enumeration.py, and it reaches r = 7.

Setting (as in cell_enumeration.py). For a tail t = (g_1, ..., g_r) of positive integers put E_0 = 0,
E_l = g_1 + ... + g_l. The comparison forms are the vectors of E_a + E_b - E_c - E_d (a + b = c + d) as linear forms in
g, one per unordered comparison, sign-normalized (first nonzero entry positive) and listed in the order produced by
`comparison_forms` (the same order as data/cells/cells_rN.json); there are 0, 1, 3, 7, 13, 22, 34 forms for
r = 1..7. Every form has coefficient sum zero, and the list contains every difference g_i - g_j.

Two-cell reduction. A real cell is given by the unit set I = {i : g_i = 1} and a full sign vector sigma of the forms.
Since the forms are translation invariant and contain all g_i - g_j, sigma fixes the set J(sigma) of minimum
coordinates; the nonempty cells with sign vector sigma are exactly C(empty, sigma) (all g_i >= 2) and
C(J(sigma), sigma), and both contain integer points (translate a witness w to w - min w + 2 and w - min w + 1).
So it suffices to enumerate the realizable sign vectors of the cone {y > 0} (equivalently g_i >= 2 after scaling and
translation) and to evaluate E-flat at both translates of each witness.

Stream format (one record per line, a final JSON summary line):
  L 0 <signs> <verdict> <w_1> ... <w_r>   a realizable full sign vector (characters + - 0 for the forms in order),
                                          with an integer witness w (w_i > 0, translated by +1 to a tail with all g_i >= 2) and the
                                          verdict code (0 Hall fails, 1 K* empty, 2 E-flat holds, 3 E-flat fails);
  F 0 <prefix> Z <j>                      an infeasible sign prefix: strict row j (in the order: the r coordinate
                                          rows y_i > 0, then the forms with sign + or - in the prefix) lies in the
                                          integer span of the equality rows (the forms with sign 0);
  F 0 <prefix> G <j>:<p>/<q> ...          an infeasible sign prefix: a positive combination of strict rows, weights
                                          summing to one, lies in the span of the equality rows (Gordan).
The terminal prefixes, read as base-3 intervals (+ = 0, - = 1, 0 = 2) of the full sign tree, must partition
[0, 3^N) in order. `verify_stream` checks every record exactly (unbounded integers), the partition, the witnesses
(every equality and strict sign), recomputes the verdict of both cells of every leaf from the integer tail by
exhaustive enumeration of row subsets (`verdict_of_tail`), and compares the counts with the summary.
"""
import gzip
import json
import lzma
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm

VERDICTS = ("nohall", "ts", "ok", "FAIL")


def level_pairs(r):
    return [[(p, s - p) for p in range(max(0, s - r), s // 2 + 1)] for s in range(2 * r + 1)]


def pair_vector(pair, r):
    """Coefficients of E_p + E_q in g_1..g_r."""
    return tuple(sum(index <= p for p in pair) for index in range(1, r + 1))


def canonical_form(a, b, r):
    vec = tuple(x - y for x, y in zip(pair_vector(a, r), pair_vector(b, r)))
    first = next(v for v in vec if v)
    sign = 1 if first > 0 else -1
    return tuple(sign * v for v in vec), sign


def comparison_forms(r):
    forms, seen = [], set()
    for group in level_pairs(r):
        for a, b in combinations(group, 2):
            vec, _ = canonical_form(a, b, r)
            if vec not in seen:
                forms.append(vec); seen.add(vec)
    return forms


def graph_verdict(attained, I, r):
    """E-flat verdict from the attaining pairs of each level, by enumerating all row subsets."""
    adj = [0] * r
    for pairs in attained:
        for p, q in pairs:
            if p < r:
                adj[p] |= 1 << q
            if q < r:
                adj[q] |= 1 << p
    unions = [0] * (1 << r); balanced = 0
    for mask in range(1, 1 << r):
        bit = mask & -mask; row = bit.bit_length() - 1
        unions[mask] = unions[mask ^ bit] | adj[row]
        size = unions[mask].bit_count()
        if size < mask.bit_count():
            return "nohall"
        if size == mask.bit_count():
            balanced |= mask
    if balanced == 0:
        return "ts"
    A = unions[balanced]
    assert A.bit_count() == balanced.bit_count()
    B = sum(1 << d for d in I if balanced & (1 << (d - 1)))
    for pairs in attained:
        if len(pairs) == 1:
            c, d = pairs[0]
            if ((A >> c) & 1 and not (B >> d) & 1) or ((A >> d) & 1 and not (B >> c) & 1):
                return "ok"
    return "FAIL"


def verdict_of_tail(t):
    r = len(t); ends = [0]
    for g in t:
        ends.append(ends[-1] + g)
    attained = []
    for group in level_pairs(r):
        mx = max(ends[p] + ends[q] for p, q in group)
        attained.append([(p, q) for p, q in group if ends[p] + ends[q] == mx])
    return graph_verdict(attained, [i + 1 for i, g in enumerate(t) if g == 1], r)


def _echelon(rows):
    """Fraction-free integer row echelon basis of the rational row span: tuple of (pivot, primitive row)."""
    red = []
    for src in rows:
        row = [F(x) for x in src]
        for piv, ex in red:
            f = row[piv]
            if f:
                row = [x - f * y for x, y in zip(row, ex)]
        piv = next((j for j, v in enumerate(row) if v), None)
        if piv is None:
            continue
        row = [v / row[piv] for v in row]
        upd = []
        for p, ex in red:
            f = ex[piv]
            if f:
                ex = [x - f * y for x, y in zip(ex, row)]
            upd.append((p, ex))
        red = sorted(upd + [(piv, row)], key=lambda it: it[0])
    out = []
    for piv, row in red:
        sc = lcm(*(v.denominator for v in row))
        vals = [int(v * sc) for v in row]
        g = gcd(*vals)
        out.append((piv, tuple(v // g for v in vals)))
    return tuple(out)


def in_span(row, ech):
    rem = list(row)
    for piv, ex in ech:
        f = rem[piv]
        if f:
            s = ex[piv]
            rem = [s * x - f * y for x, y in zip(rem, ex)]
            g = gcd(*rem)
            if g > 1:
                rem = [v // g for v in rem]
    return not any(rem)


def open_stream(path):
    p = str(path)
    return lzma.open(p, "rt") if p.endswith(".xz") else gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def verify_stream(path, r, deposited=None, second_verdict=None, second_every=1, progress=None, collect=None):
    """Check a complete stream for depth r. deposited: the cell list of data/cells/cells_rN.json (then the expanded
    cells must equal it exactly). second_verdict(tail) -> verdict: an independent verdict function applied to every
    second_every-th leaf (both cells). collect: a dict that receives sign string -> (verdict of C(empty, sigma),
    J(sigma), verdict of C(J(sigma), sigma)) for every leaf. Raises AssertionError on any failure; returns the
    counts."""
    forms = comparison_forms(r)
    for i in range(r):
        for j in range(i + 1, r):
            assert tuple(int(k == i) - int(k == j) for k in range(r)) in forms, "missing coordinate difference"
    nf = len(forms); neg = [tuple(-v for v in f) for f in forms]
    widths = [3 ** i for i in range(nf + 1)]
    trans = str.maketrans({"+": "0", "-": "1", "0": "2"})
    units = [tuple(int(i == j) for i in range(r)) for j in range(r)]

    @lru_cache(maxsize=32768)
    def eq_basis(zeros):
        return _echelon([forms[j] for j in zeros])

    expected = None
    if deposited is not None:
        assert [list(f) for f in forms] == deposited["forms"], "form order differs from the deposited cells"
        expected = {(tuple(c["unit_steps"]), tuple(c["signs"])): c["verdict"] for c in deposited["cells"]}
    found = set(); cursor = 0; counts = Counter(); cells = Counter(); kinds = Counter()
    leaves = rejected = records = second_checked = 0; summary = None
    with open_stream(path) as src:
        for line in src:
            line = line.strip()
            if not line:
                continue
            if line.startswith("{"):
                assert summary is None, "two summaries"
                summary = json.loads(line); continue
            assert summary is None, "record after the summary"
            fld = line.split(); kind = fld[0]
            assert kind in ("F", "L") and fld[1] == "0", "unknown record"
            prefix = "" if fld[2] == "_" else fld[2]
            assert len(prefix) <= nf and not prefix.strip("+-0"), "bad sign prefix"
            rank = int(prefix.translate(trans), 3) if prefix else 0
            width = widths[nf - len(prefix)]
            assert cursor == rank * width, "gap or overlap in the base-3 partition"
            cursor += width
            zeros = tuple(j for j, ch in enumerate(prefix) if ch == "0")
            strict = list(units) + [forms[j] if ch == "+" else neg[j] for j, ch in enumerate(prefix) if ch != "0"]
            if kind == "L":
                assert len(prefix) == nf, "leaf without a full sign vector"
                code = int(fld[3]); w = [int(x) for x in fld[4:]]
                assert 0 <= code < 4 and len(w) == r
                assert all(sum(a * b for a, b in zip(forms[j], w)) == 0 for j in zeros), "witness equality fails"
                assert all(sum(a * b for a, b in zip(row, w)) > 0 for row in strict), "witness sign fails"
                tail = tuple(x + 1 for x in w)
                v1 = verdict_of_tail(tail)
                assert v1 == VERDICTS[code], "verdict differs from the stream"
                mn = min(w)
                tail2 = tuple(x - mn + 1 for x in w)
                v2 = verdict_of_tail(tail2)
                assert v1 != "FAIL" and v2 != "FAIL", "E-flat fails"
                leaves += 1; counts[v1] += 1; cells[v1] += 1; cells[v2] += 1
                if collect is not None:
                    collect[prefix] = (v1, tuple(j + 1 for j, x in enumerate(w) if x == mn), v2)
                if second_verdict is not None and leaves % second_every == 0:
                    assert second_verdict(tail) == v1 and second_verdict(tail2) == v2, "second verdict differs"
                    second_checked += 1
                if expected is not None:
                    sig = tuple({"+": 1, "-": -1, "0": 0}[ch] for ch in prefix)
                    I2 = tuple(j + 1 for j, x in enumerate(w) if x == mn)
                    for key, v in (((), sig), v1), ((I2, sig), v2):
                        assert expected.get(key) == v and key not in found, "cell differs from the deposited list"
                        found.add(key)
            else:
                rejected += 1; ct = fld[3]; kinds[ct] += 1
                ech = eq_basis(zeros)
                if ct == "Z":
                    assert len(fld) == 5
                    j = int(fld[4]); assert 0 <= j < len(strict)
                    assert in_span(strict[j], ech), "zero-row certificate fails"
                elif ct == "G":
                    wts = {}
                    for tok in fld[4:]:
                        i, ratio = tok.split(":", 1); j = int(i); p, q = map(int, ratio.split("/"))
                        assert 0 <= j < len(strict) and j not in wts and p > 0 and q > 0, "bad multiplier"
                        wts[j] = (p, q)
                    den = lcm(*(q for _, q in wts.values()))
                    iw = {j: p * (den // q) for j, (p, q) in wts.items()}
                    assert sum(iw.values()) == den, "multipliers do not sum to one"
                    comb = [sum(wv * strict[j][i] for j, wv in iw.items()) for i in range(r)]
                    assert in_span(comb, ech), "Gordan certificate fails"
                else:
                    raise AssertionError("unknown certificate type")
            records += 1
            if progress and records % progress == 0:
                print(f"{records} records, {leaves} sign vectors", flush=True)
    assert summary is not None, "missing summary (incomplete stream)"
    assert cursor == widths[nf], "the terminal intervals do not cover the sign tree"
    assert summary["r"] == r and summary["forms"] == nf and (summary["unitmask_start"], summary["unitmask_stop"]) == (0, 1)
    assert summary["cells"] == leaves and summary["rejected_branches"] == rejected
    for lab, key in (("nohall", "hall_fails"), ("ts", "kstar_empty"), ("ok", "eflat_holds"), ("FAIL", "eflat_fails")):
        assert summary[key] == counts[lab], "summary counts differ"
    if expected is not None:
        assert found == set(expected), "the deposited cell list has cells not in the stream"
    return {"r": r, "forms": nf, "records": records, "sign_vectors": leaves, "rejected_prefixes": rejected,
            "rejection_kinds": dict(kinds), "cells": 2 * leaves, "cell_verdicts": dict(cells),
            "second_verdict_leaves": second_checked}
