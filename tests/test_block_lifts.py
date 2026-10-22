"""Block-lift classification of the cores (deposited data), recomputed for small k."""
import os

import pytest

from common import DATA, read_json
from dsxray.block_lifts import block_lifts


@pytest.mark.parametrize("k", [5, 6, 7, 8, 9])
def test_block_lifts_match_deposited(k):
    dep = read_json(os.path.join(DATA, "cores", f"block_lifts_k{k:02d}.json"))
    cores = read_json(os.path.join(DATA, "cores", f"cores_k{k:02d}.json"))["cores"]
    assert [r["pattern"] for r in dep["cores"]] == cores
    for r in dep["cores"]:
        assert block_lifts(tuple(r["pattern"])) == r["block_lifts"] and r["primitive"] == (not r["block_lifts"])
