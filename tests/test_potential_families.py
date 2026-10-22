"""Facet fibers, the quadratic family, refinement descent, two-point lifts, padded corners, binary cancellation and
the staircase subclasses."""
from dsxray import exact_graphs as GR
from dsxray import facet_criterion as FC
from dsxray import potential_families as PF
from dsxray import staircase_subclasses as SS
from dsxray.facet_test import classify, has_total_support


def test_exact_graphs_against_facet_criterion():
    for u in ([0, 2, 1, 3, 1, 0, 2], [0, 0, 1, 1, 0, 2, 2, 1], PF.chi(1)):
        k = len(u); _, E = GR.tight_graph(u); A = GR.allowed_edges(k, E)
        E2, A2, _ = FC.allowed_graph(u)
        assert E == set(E2) and (A is None) == (A2 is None) and (A is None or A == set(A2))


def test_cubical_fiber_m1():
    s = PF.cubical_fiber(1)
    assert (s["k"], s["lost_edges"], s["zdim_allowed"], s["bip"], s["face_dim"]) == (19, 19, 4, 1, 34)
    assert PF.zero_set_certificate()


def test_quadratic_family_small():
    for b in range(2, 5):
        for a in range(1, b + 4):
            zd, H, k = PF.quadratic(a, b, forced_edges=True)
            assert (zd == 3) == (a <= b + 1) and H == (b - 1) * a + 1
            assert classify(PF.quadratic_potential(a, b)[1])["facet"] == (a <= b + 1)


def test_descent_and_counterexample():
    for d in [(0, 1, 2), (1, 0, 2, 1), (0, 2, 1, 0)]:
        u = GR.from_increments(d)
        for m in (2, 3):
            PF.descent(u, m)
    u = GR.from_increments((0, 1, 1, 2, 0)); up = PF.spr(u, 2)
    v, E = GR.tight_graph(u)
    assert GR.zdim(len(u), E) == 4 and GR.zdim(len(up), GR.tight_graph(up)[1]) == 3
    assert not PF.lifting_condition(u, v, E) and not has_total_support(up)


def test_padding_ray_only():
    rays = [PF.padding_ray_only(w) for w in [(0, 1, 1), (1, 0, 2, 1), (0, 2, 1, 1, 0), (2, 0, 1)]]
    assert any(r[0] for r in rays)


def test_two_point_lifts_small_parent():
    w = GR.from_increments((1, 0, 0, 1))
    kids, tested = PF.two_point_lifts(w, grid=True)
    assert kids and tested > 0
    for c in kids:
        assert classify(GR.from_increments(c))["facet"]


def test_padded_corners_k20():
    words = PF.padded_corner_words(20)
    assert len(words) == 2 * 2 and len({w for *_, w in words}) == 4
    for *_, w in words:
        assert GR.zdim(20, GR.tight_graph(GR.from_increments(w))[1]) == 3


def test_binary_cancellation():
    assert PF.binary_cancellation((1, 0, 1, 1, 0, 0)) == (3, True)
    zd, cert = PF.binary_cancellation((0, 1, 0, 1))
    assert cert == (zd == 3)


def test_staircase_subclasses():
    assert SS.check_unimodal((1, 2, 3, 3, 1)) in ("Hall_fails", "strict_Hall", "endpoint_certificate")
    assert SS.no_additive_ties((4, 2, 8, 16))
    assert SS.rigid_first_row((4, 2, 8, 16), check_subsets=True) != "first_row_no_rigid_neighbor"
    assert SS.selector_counterexample() and SS.outside_range_example()
    u, want = SS.EXAMPLES["facet_without_total_support_k19"]
    assert SS.example_diagnostics(u) == want
    assert SS.bipartite_component_check([1, 2, 1, 2, 1]) in (None, "total_support", "bipartite_component")
