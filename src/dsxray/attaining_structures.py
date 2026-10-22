"""Attaining structures of sequences E_0..E_r and statement (S1); an independent checker for the catalogue.

For a sequence E (integers), M(s) = max_{p+q=s} E_p + E_q, and the attaining structure is the tuple, over the levels
s = 0..2r-1, of the sorted pairs (p, q), p <= q, p + q = s, with E_p + E_q = M(s). The classes and (S1) are as in
tools/c/attaining_structures.c and boundary_classes.py; (S1) depends on the structure only.

`check_catalogue` rechecks a deposited catalogue record by record, in exact integer arithmetic and without the C
code: the witness E realizes exactly the recorded structure (and is normalized, E_0 = E_1 = 0), the structures are
pairwise distinct, and (S1) holds on each. `brute_force_structures` lists the structures of all sequences with
small integer values, which must be a subset of a complete catalogue.
"""
from itertools import product


def pairs_at(s, r):
    return [(a, s - a) for a in range(max(0, s - r), s // 2 + 1)]


def structure_of(E):
    r = len(E) - 1; out = []
    for s in range(2 * r):
        P = pairs_at(s, r); M = max(E[a] + E[b] for a, b in P)
        out.append(tuple(p for p in P if E[p[0]] + E[p[1]] == M))
    return tuple(out)


def partners(struct, r):
    deg = [0] * (r + 1)
    for level in struct:
        for a, b in level:
            deg[a] += 1
            if b != a:
                deg[b] += 1
    return deg


def n_classes(struct, r):
    par = list(range(r))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for level in struct:
        idx = sorted({x for p in level for x in p if x < r})
        for x in idx[1:]:
            par[f(x)] = f(idx[0])
    return len({f(p) for p in range(r)})


def s1_holds(struct, r):
    deg = partners(struct, r)
    return n_classes(struct, r) <= 1 or any(deg[q] <= 1 for q in range(r))


def parse_c_line(line):
    """Parse 'AS 0:p,q;|1:...|... W e0 ... er' from tools/c/attaining_structures.c."""
    body, w = line[3:].split(" W ")
    levels = []
    for chunk in body.strip().strip("|").split("|"):
        s, rest = chunk.split(":")
        levels.append(tuple(tuple(map(int, p.split(","))) for p in rest.split(";") if p))
    return tuple(levels), [int(x) for x in w.split()]


def check_catalogue(records, r):
    """records: list of (structure, witness). Returns a dict of counts and booleans."""
    seen = set(); bad_witness = 0; s1_fail = 0; multi = 0; all2 = 0
    for struct, E in records:
        struct = tuple(tuple(tuple(p) for p in level) for level in struct)
        if len(E) != r + 1 or E[0] != 0 or E[1] != 0 or structure_of(E) != struct:
            bad_witness += 1
        seen.add(struct)
        if not s1_holds(struct, r):
            s1_fail += 1
        multi += n_classes(struct, r) > 1
        all2 += all(d >= 2 for d in partners(struct, r)[:r])
    return {"records": len(records), "distinct": len(seen), "witness_failures": bad_witness, "s1_failures": s1_fail,
            "several_classes": multi, "all_degrees_at_least_two": all2}


def brute_force_structures(r, V):
    """Structures of all integer sequences in {0..V}^(r+1)."""
    return {structure_of(E) for E in product(range(V + 1), repeat=r + 1)}
