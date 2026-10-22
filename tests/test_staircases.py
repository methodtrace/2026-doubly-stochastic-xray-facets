"""Staircases: the involution g -> g*, criterion S4, facet iff total support, rays of the fan, the three-step theorem."""
import itertools
from itertools import permutations

import pytest

from common import load_rays
from dsxray import staircases as ST
from dsxray.facet_test import has_total_support, is_facet, is_ray, v_of


@pytest.mark.parametrize("k", range(3, 10))
def test_star_is_an_involution_preserving_facets(k):
    for g in ST.staircases(k):
        assert ST.star(ST.star(g)) == g
        assert is_facet(ST.u_of(g)) == is_facet(ST.u_of(ST.star(g)))


@pytest.mark.parametrize("k", range(3, 11))
def test_S4_is_sufficient(k):
    for g in ST.staircases(k):
        if ST.S4(g):
            assert is_facet(ST.u_of(g)), g


@pytest.mark.parametrize("k", range(3, 10))
def test_facet_iff_total_support(k):
    for g in ST.staircases(k):
        assert is_facet(ST.u_of(g)) == ST.total_support(g), g


@pytest.mark.parametrize("k", [5, 6, 7, 8])
def test_fan_rays_ray_plus_total_support(k):
    for u in load_rays(k):
        u = list(u)
        assert is_ray(u) and is_facet(u) == has_total_support(u)


def test_three_step_theorem_proof_steps():
    out = ST.three_step_proof_checks(5, 14)
    for key, (n, ok) in out.items():
        assert n > 0 and ok == n, key


def test_three_step_compositions_generated_directly():
    for k in range(5, 14):
        assert sorted(g for g in ST.comps(k) if len(g) == 4) == \
            sorted((x, y - x, z - y, k - z) for x, y, z in itertools.combinations(range(1, k), 3))


def _xrays(g):
    u = ST.u_of(g); v = v_of(u); k = len(u); tot = 2 * sum(u); S = set()
    for p in permutations(range(k)):
        if sum(v[i + p[i]] for i in range(k)) == tot:
            x = [0] * (2 * k - 1)
            for i in range(k):
                x[i + p[i]] += 1
            S.add(tuple(x))
    return S


def test_star_gives_the_mirror_face():
    g = (2, 1, 1, 3); h = ST.star(g)
    assert h == (1, 1, 4, 1)
    A, B = _xrays(g), _xrays(h)
    assert A != B and {x[::-1] for x in A} == B


def test_nonincreasing_tails_small():
    for k in range(4, 11):
        for g in ST.nonincreasing_tails(k):
            assert is_facet(ST.u_of(g)) and ST.total_support(g), g
