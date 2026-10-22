"""Both vertex walks reproduce the facet walk (they share no code with it)."""
import pytest

from common import load_rays
from dsxray import vertex_walk as VW

VERTICES = {5: 22, 6: 100, 7: 484, 8: 2689}


@pytest.mark.parametrize("k", range(5, 9))
def test_exact_vertex_walk(k):
    nv, rays = VW.walk(k)
    assert nv == VERTICES[k] and rays == set(load_rays(k))


@pytest.mark.parametrize("k", range(5, 9))
def test_cdd_vertex_walk(k):
    pytest.importorskip("cdd")
    from dsxray import vertex_walk_cdd as VC
    nv, rays = VC.walk(k)
    assert nv == VERTICES[k] and rays == set(load_rays(k))


def test_cdd_facets_k6():
    pytest.importorskip("cdd")
    from dsxray.small_k_fan import cdd_facets_verified, exact_facets
    assert cdd_facets_verified(6) == exact_facets(6)[0]
