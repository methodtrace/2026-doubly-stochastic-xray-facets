"""Padding theorem on deposited facets, cores, and the families A, W, V."""
import os

import pytest

from common import DATA, load_rays, read_json
from dsxray import nonstep_families as NF
from dsxray import padding as PD
from dsxray.facet_test import classify, tight_edges


def facets(k):
    return {PD.pattern(u) for u in load_rays(k) if classify(list(u))["facet"]}


@pytest.mark.parametrize("k", [6, 7])
def test_padding_maps_facets_to_facets(k):
    F, F1 = facets(k), facets(k + 1)
    assert all(PD.pad_left(p) in F1 and PD.pad_right(p) in F1 for p in F)


def test_padding_hypotheses_and_structure_k8():
    for p in facets(8):
        ts, conn = PD.connected_with_total_support(PD.from_pattern(p))
        assert ts and conn and PD.P_L(p) and PD.P_R(p) and PD.padding_structure_ok(p)


def test_cores_up_to_8():
    F = {3: PD.small_facet_patterns(3), 4: PD.small_facet_patterns(4)}
    F.update({k: facets(k) for k in range(5, 9)})
    table, cores = PD.core_table(F, 3, 8)
    for k in range(3, 9):
        assert cores[k] == sorted(map(tuple, read_json(os.path.join(DATA, "cores", f"cores_k{k:02d}.json"))["cores"]))
        assert len(F[k]) == sum(len(cores[j]) * (k - j + 1) for j in range(3, k + 1))
    assert (table[8]["cores"], table[8]["nonstep_cores"]) == (62, 36)


def test_mirror_padding_relation():
    p = (1, 0, 2, 0, 1)
    assert PD.pad_right(p) == PD.mirror(PD.pad_left(PD.mirror(p)))


@pytest.mark.parametrize("k", range(9, 17))
def test_family_proof_objects(k):
    assert tight_edges(NF.W_u(k)) == NF.W_tight_graph(k)
    assert NF.permutation_with_universal_vertex(NF.W_tight_graph(k), *NF.W_permutation(k))
    assert tight_edges(NF.Vplus_u(k)) == NF.Vplus_tight_graph(k)
    assert NF.permutation_with_universal_vertex(NF.Vplus_tight_graph(k), *NF.Vplus_permutation(k))
    assert tight_edges(NF.Vminus_u(k)) == NF.Vminus_tight_graph(k)
    assert NF.permutation_with_universal_vertex(NF.Vminus_tight_graph(k), *NF.Vminus_permutation(k))


@pytest.mark.parametrize("k", range(6, 13))
def test_family_membership(k):
    for a in range(2, k - 3):
        u = NF.family_A_u(a, k - 3 - a)
        assert NF.check_member(u) == (True, True, True)
        assert tight_edges(u) == NF.family_A_tight_graph(a, k - 3 - a)
    assert NF.check_member(NF.W_u(k))[0] == (k == 7 or k >= 9)
    assert NF.check_member(NF.Vplus_u(k))[0] == (k >= 7)
    assert NF.check_member(NF.Vminus_u(k))[0] == (k == 7 or k >= 9)
