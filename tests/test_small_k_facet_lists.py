"""Exact facet lists for k <= 8 from all permutation profiles, line enumeration of the fan, and agreement with the
certified test on the deposited rays."""
import pytest

from common import load_rays
from dsxray.facet_test import classify, classify_coordinate
from dsxray.small_k_fan import coordinate_inequality, exact_facets, normalize_ray, profiles, ray_inequality

EXPECTED = {2: 2, 3: 6, 4: 13, 5: 23, 6: 45, 7: 85, 8: 187}


@pytest.mark.parametrize("k", range(2, 9))
def test_exact_facet_counts(k):
    H, rays = exact_facets(k)
    assert len(H) == EXPECTED[k]
    if k >= 5:
        assert {normalize_ray(u) for u in rays} == set(load_rays(k))
        mine = {ray_inequality(u) for u in load_rays(k) if classify(list(u))["facet"]}
        mine |= {coordinate_inequality(k, t) for t in range(2 * k - 1) if classify_coordinate(k, t)["facet"]}
        assert mine == H


def test_profile_count_k8():
    assert len(profiles(8)) == 5270
