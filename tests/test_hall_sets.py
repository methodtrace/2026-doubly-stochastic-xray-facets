"""thm:H: matching criterion, brute force over sets of starts, total support, and the proof inequalities."""
import pytest

from dsxray import hall_sets as HS
from dsxray.staircases import staircases


@pytest.mark.parametrize("k", range(4, 11))
def test_theorem_and_brute_force(k):
    c = HS.theorem_counts(k, brute=True)
    assert c["exceptions"] == 0 and c["brute_force_disagreements"] == 0 and c["neighbors_not_ends"] == 0


def test_proof_inequalities():
    n = sum(HS.proof_step_checks(g) for k in range(3, 10) for g in staircases(k))
    assert n == 6026 + 15249


def test_minimal_sizes_small():
    n, sizes = HS.minimal_size_counts(3, 10)
    assert n > 0 and set(sizes) <= {"1", "2", "3"}
