"""(S1), lem:R and the (H+) counterexamples (not used in the paper); agreement with the deposited exhaustive counts."""
import os

import pytest

from common import DATA, read_json
from dsxray import boundary_classes as BC

DEP = {r["k"]: r for r in read_json(os.path.join(DATA, "summary", "boundary_classes.json"))["s1_compositions"]}


@pytest.mark.parametrize("k", range(4, 15))
def test_s1_matches_deposited_counts(k):
    c = BC.s1_counts(k)
    assert c["s1_failures"] == 0 and c == DEP[k]


def test_s1_pattern_lifts_small():
    for V, r in ((2, 6), (3, 5)):
        assert BC.s1_pattern_lifts(V, r)["s1_failures"] == 0


def test_lemma_r():
    c = BC.lemma_r_counts(4, 12)
    assert c["with_solitary_index"] > 0 and c["of_which_total_support"] == 0


def test_hplus_counterexample():
    d = BC.describe((1, 1, 4, 4, 7, 4, 4, 1, 6, 5))
    assert d["k"] == 37 and d["hplus"] == "FAIL" and not d["total_support"] and d["dim_z"] == 3 and d["s1_holds"]
    assert d["classes"] == [[0, 1, 2, 3, 4, 5, 6, 7], [8]] and d["solitary"] == [7]
    assert (d["unused_tight_edges"], d["tight_edges"]) == (362, 399)
