"""Facet walk: ray counts, facet counts, consistency, parallel = serial, deposited lists, census and arithmetic bound."""
import pytest

from common import load_facet_rows, load_rays
from dsxray import fan_walk as FW
from dsxray.facet_test import classify, classify_coordinate

RAYS = {5: 16, 6: 42, 7: 110, 8: 338, 9: 1146}
FACETS = {5: 23, 6: 45, 7: 85, 8: 187, 9: 426}


@pytest.mark.parametrize("k", range(5, 10))
def test_walk(k):
    rays, st = FW.fan_rays(k)
    assert len(rays) == RAYS[k] and st["adjacency_symmetric"] and set(rays) == set(load_rays(k))
    nf = sum(classify(list(u))["facet"] for u in rays) + sum(classify_coordinate(k, t)["facet"] for t in range(2 * k - 1))
    assert nf == FACETS[k]


def test_parallel_walk_equals_serial():
    r1, _ = FW.fan_rays(8); r2, st = FW.fan_rays_parallel(8, 1)
    assert r1 == r2 and st["adjacency_symmetric"]


def test_recursive_variant_k8():
    rays, st = FW.fan_rays(8, walk_rank=12)
    assert set(rays) == set(load_rays(8)) and st["walks"] > 0


def test_resumable_walk(tmp_path):
    rays, st = FW.fan_rays_resumable(7, 1, str(tmp_path / "state.jsonl"), walk_rank=None)
    assert set(rays) == set(load_rays(7)) and st["one_sided_ridges"] == 0
    rays2, _ = FW.fan_rays_resumable(7, 1, str(tmp_path / "state.jsonl"), walk_rank=None)
    assert rays2 == rays


def test_rank_mod_p_is_exact_in_range():
    assert FW.hadamard_ok(25) and not FW.hadamard_ok(26)
    for k in (6, 9):
        norms = {sum(x * x for x in FW.dvec(a, b, k)) for S in FW.Q_summands(k) for a in S for b in S if a != b}
        assert max(norms) == 6


def test_census_k9():
    F = [r["u"] for r in load_facet_rows(9) if r["facet"]]
    c = FW.census(9, F, 17)
    assert c == {"coordinate": 17, "single": 16, "double": 56, "triple": 88, "four-step": 54, "other": 195}
