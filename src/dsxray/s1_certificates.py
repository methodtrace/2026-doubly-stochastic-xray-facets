"""Statement (S1) for all real sequences with r steps: exact verification of a certificate tree.

For a real sequence E_0..E_r put M(s) = max_{p+q=s} E_p + E_q. At each level s = 0..2r-1 the attaining set is a
nonempty subset of the pairs (p, s - p), p <= s - p (`pairs_at`); the level 2r (only the pair (r, r)) affects neither
the boundary classes nor the partner sets of q < r and is omitted. The boundary classes are the connected components
on {0..r-1} of the relation "both are entries < r of attaining pairs at a common level"; P(q) is the set of partners
of q in attaining pairs. (S1): if there are at least two classes, some q < r has |P(q)| <= 1.

Adding alpha + beta p to E_p changes no attaining set, so E_0 = E_1 = 0 may be imposed. The attaining sets realized by
some sequence form a finite tree of choices (level by level, a nonempty subset of the pairs, encoded as a mask with
radix 2^(number of pairs) - 1); the sequences realizing a prefix of choices form the relatively open cone
{E : E_0 = E_1 = 0, equalities between the chosen pairs, strict inequalities "first chosen pair - unchosen pair > 0"}.
A certificate tree lists terminal prefixes (masks 1..radix, comma separated, `_` for the empty prefix):
  C <prefix>            the boundary classes of the chosen pairs are already connected (later levels only merge);
  D <prefix> <q>        q < r has fewer than two possible partners: chosen ones plus every l with q + l at a level
                        not yet chosen;
  Z <prefix> <j>        strict row j lies in the row span of the equalities (the cone is empty);
  G <prefix> <j>:<p>/<q> ...  a positive combination of strict rows (weights summing to one) lies in that span
                        (Gordan; the cone is empty);
  X <prefix> <E_0..E_r> a counterexample leaf (a sequence realizing the full prefix without (S1)); none occur.
The first line `T {...}` gives r and the domain prefix, the last line a JSON summary with the counts. The terminal
prefixes, read as mixed-radix intervals, must partition the domain in order. `verify_tree` checks all this in exact
integer arithmetic. `probes` are sequences: for each, the terminal interval containing its full attaining structure
must be a C or D record whose claim holds for the actual sequence (an end-to-end test of the encoding).
"""
import gzip
import hashlib
import json
import lzma
from collections import Counter
from functools import lru_cache
from math import gcd, lcm

from .attaining_structures import n_classes, pairs_at, partners, structure_of


def _primitive(row):
    g = gcd(*row)
    return tuple(v // g for v in row) if g > 1 else tuple(row)


def _echelon(rows):
    out = []
    for src in rows:
        row = tuple(src)
        for piv, ex in out:
            f = row[piv]
            if f:
                row = _primitive([ex[piv] * a - f * b for a, b in zip(row, ex)])
        piv = next((j for j, v in enumerate(row) if v), None)
        if piv is None:
            continue
        out = sorted(out + [(piv, _primitive(row))], key=lambda it: it[0])
    return tuple(out)


def _in_span(src, ech):
    row = tuple(src)
    for piv, ex in ech:
        f = row[piv]
        if f:
            row = _primitive([ex[piv] * a - f * b for a, b in zip(row, ex)])
    return not any(row)


def parse_prefix(value):
    return () if value == "_" else tuple(int(x) for x in value.split(","))


def tree_shape(r):
    groups = [pairs_at(s, r) for s in range(2 * r)]
    radices = [(1 << len(g)) - 1 for g in groups]
    widths = [1] * (2 * r + 1)
    for s in range(2 * r - 1, -1, -1):
        widths[s] = radices[s] * widths[s + 1]
    return groups, radices, widths


def interval(prefix, radices, widths):
    rank = 0
    for s, mask in enumerate(prefix):
        assert 1 <= mask <= radices[s], "invalid mask"
        rank = rank * radices[s] + mask - 1
    return rank * widths[len(prefix)], (rank + 1) * widths[len(prefix)]


def structure_rank(E, r):
    """Position of the full attaining structure of E (levels 0..2r-1) in the mixed-radix tree."""
    groups, radices, widths = tree_shape(r)
    st = structure_of(E); masks = []
    for s, lev in enumerate(st):
        masks.append(sum(1 << groups[s].index(p) for p in lev))
    return interval(tuple(masks), radices, widths)[0]


def _components(levels, n):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for lev in levels:
        idx = sorted({x for p in lev for x in p if x < n})
        for x in idx[1:]:
            par[f(x)] = f(idx[0])
    return len({f(p) for p in range(n)})


def open_binary(path):
    p = str(path)
    return lzma.open(p, "rb") if p.endswith(".xz") else gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")


def verify_tree(path, r, domain_prefix=(), probes=(), progress=None):
    """Check one certificate tree (or one domain of it); path is a file name or a function returning a binary file
    object of the uncompressed stream. Raises AssertionError on failure; returns counts, the
    SHA-256 and size of the uncompressed stream, and the number of probes checked."""
    groups, radices, widths = tree_shape(r)
    start, end = interval(tuple(domain_prefix), radices, widths); cursor = start
    units = [tuple(int(i == j) for i in range(r + 1)) for j in (0, 1)]
    pr = sorted((structure_rank(E, r), tuple(E)) for E in probes)
    pr = [x for x in pr if start <= x[0] < end]; pi = 0; probed = 0

    def diff(a, b):
        row = [0] * (r + 1)
        for p in a:
            row[p] += 1
        for p in b:
            row[p] -= 1
        return tuple(row)

    @lru_cache(maxsize=16384)
    def cons(prefix):
        if not prefix:
            return tuple(units), ()
        EQ, IN = cons(prefix[:-1]); s = len(prefix) - 1; mask = prefix[-1]; g = groups[s]
        first = next(i for i in range(len(g)) if mask >> i & 1); eq = list(EQ); st = list(IN)
        for i, pair in enumerate(g):
            if i != first:
                (eq if mask >> i & 1 else st).append(diff(g[first], pair))
        return tuple(eq), tuple(st)

    @lru_cache(maxsize=32768)
    def ech(EQ):
        return _echelon(EQ)

    def chosen(prefix):
        return tuple(tuple(p for i, p in enumerate(groups[s]) if m >> i & 1) for s, m in enumerate(prefix))

    counts = Counter(); header = summary = None; h = hashlib.sha256(); nbytes = 0; records = 0
    with (path() if callable(path) else open_binary(path)) as src:
        for raw in src:
            h.update(raw); nbytes += len(raw)
            line = raw.decode().strip()
            if not line:
                continue
            if line.startswith("T "):
                assert header is None and records == 0, "misplaced header"
                header = json.loads(line[2:])
                assert header["r"] == r and tuple(header["domain_prefix"]) == tuple(domain_prefix), "wrong header"
                continue
            assert header is not None, "missing header"
            if line.startswith("{"):
                assert summary is None, "two summaries"
                summary = json.loads(line); continue
            assert summary is None, "record after the summary"
            f = line.split(); kind = f[0]
            assert kind in "CDZGX" and len(kind) == 1, "unknown record"
            prefix = parse_prefix(f[1]); lo, hi = interval(prefix, radices, widths)
            lo, hi = max(lo, start), min(hi, end)
            assert lo < hi and cursor == lo, "gap or overlap in the partition"
            cursor = hi
            if kind == "C":
                assert len(f) == 2 and _components(chosen(prefix), r) == 1, "invalid connectivity cut"
            elif kind == "D":
                assert len(f) == 3; q = int(f[2]); assert 0 <= q < r
                part = {b if a == q else a for lev in chosen(prefix) for a, b in lev if q in (a, b)}
                fut = {l for l in range(r + 1) if q + l >= len(prefix)}
                assert len(part | fut) < 2, "invalid solitary-index cut"
            else:
                EQ, IN = cons(prefix); B = ech(EQ)
                if kind == "Z":
                    assert len(f) == 3; j = int(f[2]); assert 0 <= j < len(IN)
                    assert _in_span(IN[j], B), "zero-row certificate fails"
                elif kind == "G":
                    w = {}
                    for tok in f[2:]:
                        i, ratio = tok.split(":", 1); j = int(i); p, q = map(int, ratio.split("/"))
                        assert 0 <= j < len(IN) and j not in w and p > 0 and q > 0, "bad multiplier"
                        w[j] = (p, q)
                    assert w
                    den = lcm(*(q for _, q in w.values()))
                    iw = {j: p * (den // q) for j, (p, q) in w.items()}
                    assert sum(iw.values()) == den, "multipliers do not sum to one"
                    comb = [sum(x * IN[j][i] for j, x in iw.items()) for i in range(r + 1)]
                    assert _in_span(comb, B), "Gordan certificate fails"
                else:
                    raise AssertionError("counterexample leaf in the tree")
            while pi < len(pr) and pr[pi][0] < hi:
                rank, E = pr[pi]
                assert lo <= rank and kind in "CD", "a realized structure lies in an infeasible branch"
                st = structure_of(list(E))
                if kind == "C":
                    assert n_classes(st, r) == 1
                else:
                    assert partners(st, r)[int(f[2])] <= 1
                pi += 1; probed += 1
            counts[kind] += 1; records += 1
            if progress and records % progress == 0:
                print(f"{records} terminal records", flush=True)
    assert summary is not None, "missing summary (incomplete stream)"
    assert cursor == end, "the terminals do not cover the domain"
    assert summary["mode"] == "s1-proof" and summary["r"] == r
    assert summary["connectivity_pruned"] == counts["C"] and summary["degree_pruned"] == counts["D"]
    assert summary["rejected_branches"] == counts["Z"] + counts["G"] and summary["infeasibility_certificates"] == counts["G"]
    assert summary["leaves"] == counts["X"] == 0 and summary["s1_failures"] == 0
    return {"r": r, "domain_prefix": list(domain_prefix), "domain_interval": [start, end], "terminal_records": records,
            "certificate_types": dict(counts), "uncompressed_sha256": h.hexdigest(), "uncompressed_bytes": nbytes,
            "probes_checked": probed}
