"""Facet criterion (thm:F), empty anti-diagonals (thm:empty), symmetric dual (thm:symdual), ex:noTS."""
import os
import random

import pytest

from common import DATA, load_rays, read_json
from dsxray.facet_criterion import (allowed_graph, bip, check_symmetric_dual, components, criterion, rank_int,
                                    symmetric_dual, y_dim)
from dsxray.facet_test import classify, rank_exact

POT = {p["name"]: p["u"] for p in read_json(os.path.join(DATA, "counterexamples", "potentials.json"))["potentials"]}


def test_rank_int_matches_fraction_rank():
    rng = random.Random(4)
    for _ in range(60):
        rows = [[rng.randint(-2, 2) for _ in range(rng.randint(1, 7))] for _ in range(rng.randint(1, 6))]
        n = max(len(r) for r in rows); rows = [r + [0] * (n - len(r)) for r in rows]
        assert rank_int(rows) == rank_exact(rows)


@pytest.mark.parametrize("k", [6, 7])
def test_criterion_equals_facet_test_on_rays(k):
    for u in load_rays(k):
        u = list(u); c = criterion(u, face_dim=True)
        assert c["criterion_facet"] == classify(u)["facet"]
        if c["tight_permutation"]:
            E, A, _ = allowed_graph(u)
            assert y_dim(k, A) == c["bip"] == bip(k, A)
            if c["cover"]:
                assert c["face_dim"] == c["formula_dim"]
            else:
                assert c["face_dim"] < 2 * k - 4


def test_example_5_1_and_remark_5_2():
    chi, chp = POT["example_5_1_chi"], POT["remark_5_2_chi_prime"]
    for u in (chi, chp):
        c = criterion(u, face_dim=True)
        assert c["zdim"] == 3 and not c["total_support"] and c["cover"] and c["bip"] == 1 and c["zdim_star"] == 4
        assert c["face_dim"] == 34 and classify(u)["facet"]
    E, A, _ = allowed_graph(chi)
    assert ([11, 16], True) in components(19, A) and len(E) - len(A) == 19


def test_symmetric_dual_examples():
    for name in ("example_5_1_chi", "remark_5_2_chi_prime", "k17_cover_without_total_support"):
        u = POT[name]; w = symmetric_dual(u)
        assert w is not None and w != u
        assert all(check_symmetric_dual(u, w).values())


def test_k17_bipartite_component_with_total_support():
    W = POT["k17_total_support_bipartite_component"]; c = criterion(W)
    E, A, _ = allowed_graph(W)
    assert c["total_support"] and c["zdim"] == 6 and ([8, 10], True) in components(17, E)
