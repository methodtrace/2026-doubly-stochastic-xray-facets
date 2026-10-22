"""E-flat: cell enumeration for small r, integer tails, and consistency with the tight graph."""
import os

import pytest

from common import DATA, read_json
from dsxray import cell_enumeration as CE
from dsxray.staircases import staircases

EXPECTED = {1: (2, 0, 2, 0), 2: (6, 0, 6, 0), 3: (26, 0, 20, 6), 4: (198, 12, 90, 96)}


@pytest.mark.parametrize("r", range(1, 5))
def test_cells(r):
    C = CE.CellEnumeration(r); C.run(); s = C.summary()
    assert (s["cells"], s["hall_fails"], s["kstar_empty"], s["eflat_holds"]) == EXPECTED[r] and s["eflat_fails"] == 0
    n, realized, miss = C.safeguard(14)
    assert miss == 0
    dep = read_json(os.path.join(DATA, "cells", f"cells_r{r}.json"))
    assert dep == C.export()


def test_eflat_tails():
    dep = {r["n"]: r["counts"] for r in read_json(os.path.join(DATA, "summary", "eflat.json"))["eflat_tails"]}
    for n in range(1, 13):
        assert CE.eflat_counts(n)["counts"] == dep[n]


def test_feasibility_small():
    from fractions import Fraction as F
    assert CE.feasible([((F(1),), F(-1), '>'), ((F(-1),), F(2), '>')], 1)       # 1 < x < 2
    assert not CE.feasible([((F(1),), F(-2), '>'), ((F(-1),), F(2), '>=')], 1)  # x > 2 and x <= 2


def test_consistency_with_tight_graph():
    for k in range(4, 11):
        for g in staircases(k):
            CE.largest_closed_start_set_matches(g)
