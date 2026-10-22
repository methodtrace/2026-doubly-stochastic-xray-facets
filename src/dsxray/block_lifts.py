"""Block faces of hat-P(k) and the block-lift classification of cores (deposited data; not cited in the paper).

For 1 <= b <= k - 3 the inequality sum_{t < b} x(t) <= b is valid on hat-P(k), and equality forces the matrix entries
Z(i, b-1-i) = 1 for i < b: the face F_L(b) is the set of X-rays of permutations that reverse {0..b-1}, affinely
isomorphic to hat-P(k - b). F_R(b) is its mirror image (the last b indices reversed). A facet ray u is a block lift
on side s at b if its facet meets F_s(b) in a facet of F_s(b), that is, if the permutations tight for u that reverse
the block span a face of dimension 2(k - b) - 4. A core with no block lift is called primitive.

`block_lifts(p)` returns the list of (side, b) for the normalized pattern p; face dimensions are exact (allowed edges
of the restricted tight graph, then rank [M; X] - rank M with exact integer ranks).
"""
from .facet_criterion import face_dim_exact
from .facet_test import allowed_edges, perfect_matching, tight_edges
from .padding import from_pattern


def block_face_edges(E, k, b, side):
    K = k - 1; out = set()
    for (i, j) in E:
        if side == "L":
            ib, jb = i < b, j < b
            if (ib and jb and i + j == b - 1) or (not ib and not jb):
                out.add((i, j))
        else:
            ib, jb = i > K - b, j > K - b
            if (ib and jb and i + j == 2 * K - (b - 1)) or (not ib and not jb):
                out.add((i, j))
    return out


def face_dim_of_edges(k, edges):
    if not edges:
        return -1
    m = perfect_matching(k, edges)
    if m is None:
        return -1
    return face_dim_exact(k, allowed_edges(k, edges, m))


def block_lifts(p):
    u = from_pattern(p); k = len(u); E = tight_edges(u); out = []
    for b in range(1, k - 2):
        for side in "LR":
            if face_dim_of_edges(k, block_face_edges(E, k, b, side)) == 2 * (k - b) - 4:
                out.append([side, b])
    return out
