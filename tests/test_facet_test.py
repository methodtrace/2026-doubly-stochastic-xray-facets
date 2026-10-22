"""The certified facet test agrees with two independent tests (allowed edges by deletion; enumeration of all
permutation profiles) and reproduces the known families, examples and coordinate facets."""
import itertools

import pytest

from dsxray.facet_test import (classify, classify_coordinate, has_total_support, is_facet, reference_face_dim,
                               reference_is_facet, matchable_edges, rank_exact, tight_edges, face_rows, v_of)
from dsxray.small_k_fan import profile_is_facet, profiles


def u_of(J, k):
    return [sum(int(i >= j) for j in J) for i in range(k)]


FAMILIES = {
    "single": lambda k: [u_of((j,), k) for j in range(1, k)],
    "double": lambda k: [u_of((a, b), k) for a in range(1, k) for b in range(a + 1, k)],
    "triple": lambda k: [u_of(J, k) for J in itertools.combinations(range(1, k), 3)],
    "hinge": lambda k: [[max(0, i - j) for i in range(k)] for j in range(0, k - 1)],
    "floor": lambda k: [[i // d for i in range(k)] for d in range(2, k)],
}


@pytest.mark.parametrize("k", [3, 4, 5, 6, 7])
def test_three_tests_agree(k):
    P = profiles(k)
    for name, gen in FAMILIES.items():
        for u in gen(k):
            if len(set(u[i + 1] - u[i] for i in range(k - 1))) <= 1:
                continue                                  # affine potentials are not considered
            for sg in (1, -1):
                f = is_facet(u, sg)
                assert f == reference_is_facet(u, sg) == profile_is_facet(u, sg, P), (name, u, sg)


@pytest.mark.parametrize("k,fails", [(6, 0), (7, 6), (8, 12), (9, 24)])
def test_triple_step_failure_counts(k, fails):
    assert sum(not is_facet(u, sg) for u in FAMILIES["triple"](k) for sg in (1, -1)) == fails


@pytest.mark.parametrize("k", [10, 11, 12])
def test_single_and_double_steps_are_facets(k):
    for name in ("single", "double"):
        for u in FAMILIES[name](k):
            assert is_facet(u, 1) and is_facet(u, -1)


@pytest.mark.parametrize("k", range(4, 10))
def test_consecutive_steps_are_facets(k):
    for r in range(1, k - 1):
        for a in range(1, k - r + 1):
            for sg in (1, -1):
                assert is_facet(u_of(tuple(range(a, a + r)), k), sg)


@pytest.mark.parametrize("k", range(5, 10))
def test_equal_gap_and_top_step_triples(k):
    for J in itertools.combinations(range(1, k), 3):
        if J[1] - J[0] == J[2] - J[1]:
            assert is_facet(u_of(J, k), 1) and is_facet(u_of(J, k), -1)
        if J[2] == k - 1:
            assert is_facet(u_of(J, k), 1)


@pytest.mark.parametrize("k", range(7, 11))
def test_triple_1_2_k_minus_2_depends_on_sign(k):
    assert not is_facet(u_of((1, 2, k - 2), k), 1) and is_facet(u_of((1, 2, k - 2), k), -1)


@pytest.mark.parametrize("k", range(4, 8))
def test_sign_identity(k):
    """(J, -) and (J^c, +) differ by an affine function, so they give the same inequality."""
    for r in range(1, k - 1):
        for J in itertools.combinations(range(1, k), r):
            Jc = tuple(j for j in range(1, k) if j not in J)
            d = [-x - y for x, y in zip(u_of(J, k), u_of(Jc, k))]
            assert all(d[i + 2] - 2 * d[i + 1] + d[i] == 0 for i in range(k - 2))


def test_raw_tight_graph_would_give_a_false_facet():
    """For (1, 2, 5) at k = 7 the raw tight graph has the facet dimension; the allowed graph has dimension 1."""
    k = 7; u = u_of((1, 2, 5), k); E = sorted(tight_edges(u)); Ep = sorted(matchable_edges(k, set(E)))

    def dim(ed):
        M, A = face_rows(k, ed)
        return rank_exact(M + A) - rank_exact(M)
    assert dim(E) == 2 * k - 4 and dim(Ep) == 1 and reference_face_dim(u, exact=True)[0] == 1
    assert not classify(u)["facet"]


def test_empty_face():
    c = classify(u_of((1, 2, 5), 10))
    assert c["cert"] == "empty" and not c["facet"]


def test_coordinate_facets():
    assert [classify_coordinate(3, t)["facet"] for t in range(5)] == [True, True, False, True, True]
    for k in range(4, 10):
        assert all(classify_coordinate(k, t)["facet"] for t in range(2 * k - 1))


def test_hinge_formulas_and_facets():
    for k in range(3, 14):
        for c in range(1, k - 1):
            assert v_of([max(i - c, 0) for i in range(k)]) == [max(t - 2 * c, 0) for t in range(2 * k - 1)]
            u = [-max(i - c, 0) for i in range(k)]
            fold = [(-max(t - c, 0)) if t <= k - 1 else (-max(t - k + 1 - c, 0) - (k - 1 - c)) for t in range(2 * k - 1)]
            assert v_of(u) == fold
    for k in range(4, 12):
        for c in range(1, k - 1):
            for sg in (1, -1):
                assert is_facet([max(i - c, 0) for i in range(k)], sg)


def test_total_support_helper():
    assert has_total_support(u_of((3,), 8)) and not has_total_support(u_of((1, 2, 5), 7))
