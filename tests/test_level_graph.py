"""Level structure of staircases: starts, last anti-diagonals, square relations, (H+) for small k, rigid corners."""
import pytest

from dsxray import level_graph as LG
from dsxray import staircases as ST
from dsxray.facet_test import classify

RIGID = {7: (59, 3, 0), 8: (114, 12, 0), 9: (214, 37, 3), 10: (392, 104, 14), 11: (710, 254, 58),
         12: (1270, 600, 176), 13: (2245, 1324, 525)}


def test_block_start_neighbors_and_corners():
    for k in range(3, 12):
        for g in ST.staircases(k):
            assert LG.start_neighbors_are_attaining_ends(g) and LG.last_antidiagonal_corners(g), g


def test_gstar_reflects_tight_graph():
    for k in range(3, 11):
        for g in ST.staircases(k):
            assert LG.gstar_tight_graph_reflects(g), g


def test_k22_example():
    g = (1, 1, 5, 3, 2, 2, 3, 5)
    a = LG.square_relation_analysis(g)
    assert a["total_support"] and a["ray"] and not a["boundary_connected"]
    assert LG.hplus_verdict(list(g))[0] == "connected" and classify(ST.u_of(g))["facet"]


@pytest.mark.parametrize("k", range(4, 15))
def test_hplus_small_k(k):
    assert "FAIL" not in LG.hplus_counts(k)["counts"]


def test_one_class_under_total_support():
    for k in range(4, 13):
        for g in ST.staircases(k):
            if ST.total_support(g):
                assert LG.hplus_verdict(list(g))[0] == "connected" and LG.hplus_verdict(list(ST.star(g)))[0] == "connected"


@pytest.mark.parametrize("k", sorted(RIGID))
def test_rigid_corner(k):
    C = {}; bad = 0
    for g in ST.staircases(k):
        s = LG.rigid_corner(g); C[s] = C.get(s, 0) + 1
        bad += LG.first_single_start_partner(g) is False
    assert "FAIL" not in C and bad == 0
    assert (C.get("total_support", 0), C.get("certified", 0), C.get("no_pm", 0)) == RIGID[k]
