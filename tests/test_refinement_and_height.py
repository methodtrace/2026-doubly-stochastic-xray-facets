"""Grid refinement (thm:refine; lem:sprtight, lem:blockgraph, lem:flatlevel) and the family G_{1,b} of thm:height."""
import pytest

from common import load_rays
from dsxray import nonstep_families as NF
from dsxray import padding as PD
from dsxray import refinement as RF
from dsxray.facet_criterion import graph_z_dim
from dsxray.facet_test import classify, tight_edges, v_of


def facet_patterns(k):
    return sorted({PD.pattern(u) for u in load_rays(k) if classify(list(u))["facet"]})


@pytest.mark.parametrize("m", [2, 3])
def test_refinement_k6(m):
    for p in facet_patterns(6):
        u = PD.from_pattern(p)
        assert RF.tight_graph_lemma_holds(u, m) and RF.theta_connected(u) and RF.flat_level(u) is not None
        up = RF.spr(u, m)
        assert graph_z_dim(len(up), tight_edges(up)) == 3 and classify(up)["facet"]
        assert PD.connected_with_total_support(up) == (True, True) and PD.P_L(PD.pattern(up))
        assert classify(RF.ref(u, m))["facet"]


def test_ref_is_padding_of_spr():
    u = PD.from_pattern((1, 0, 2, 0, 1))
    for m in (2, 3, 4):
        q = PD.pattern(RF.spr(u, m))
        for _ in range(m - 1):
            q = PD.pad_left(q)
        assert PD.pattern(RF.ref(u, m)) == q


@pytest.mark.parametrize("b", [2, 3, 4, 5, 6])
def test_height_family(b):
    p = NF.height_family_pattern(b); u = PD.from_pattern(p); k = len(u)
    w = [u[i] - i for i in range(k)]
    assert v_of(w) == NF.height_family_v(b) and tight_edges(u) == NF.height_family_tight_graph(b)
    s = NF.height_family_permutation(b)
    assert sorted(s) == list(range(k)) and all((i, s[i]) in tight_edges(u) for i in range(k))
    assert classify(u)["facet"] and max(PD.pattern(u)) == b and PD.P_L(p) and PD.P_R(p)


def test_quadratic_height_family_small():
    for b in (2, 3):
        for a in range(1, b + 4):
            assert classify(PD.from_pattern(NF.quadratic_height_pattern(a, b)))["facet"] == (a <= b + 1)
