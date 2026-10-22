"""The exhaustive search program for failures of total support, against Python on small cases."""
import itertools
import json
import os
import shutil
import subprocess

import pytest

from common import ROOT
from dsxray.facet_criterion import criterion
from dsxray.padding import from_pattern


@pytest.fixture(scope="module")
def exe(tmp_path_factory):
    cc = shutil.which("cc") or shutil.which("gcc")
    if not cc:
        pytest.skip("no C compiler")
    out = str(tmp_path_factory.mktemp("bin") / "search_total_support")
    subprocess.run([cc, "-O2", "-o", out, os.path.join(ROOT, "tools", "c", "search_total_support.c")], check=True)
    return out


def run(exe, *args):
    lines = subprocess.run([exe] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout.splitlines()
    return json.loads(lines[-1]), lines[:-1]


def test_values_mode_against_python(exe):
    k, H = 6, 2
    summary, _ = run(exe, "V", k, H)
    tp = ts = 0; hist = {}
    for rest in itertools.product(range(H + 1), repeat=k - 1):
        c = criterion([0] + list(rest))
        if c["tight_permutation"]:
            tp += 1; ts += c["total_support"]
            if not c["total_support"]:
                hist[str(c["n_empty"])] = hist.get(str(c["n_empty"]), 0) + 1
    assert (summary["tight_permutation"], summary["total_support"]) == (tp, ts) and summary["histogram_empty"] == hist


def test_increments_mode_against_python(exe):
    k, H = 7, 2
    summary, found = run(exe, "I", k, H)
    rays = nots = 0
    for d in itertools.product(range(H + 1), repeat=k - 1):
        if min(d) != 0 or max(d) == 0 or (max(d) == 2 and all(x % 2 == 0 for x in d)):
            continue
        c = criterion(from_pattern(d))
        if c["tight_permutation"] and c["zdim"] == 3:
            rays += 1; nots += not c["total_support"]
    assert (summary["rays_with_tight_permutation"], summary["rays_without_total_support"]) == (rays, nots) and not found


def test_finds_example_5_1(exe):
    summary, lines = run(exe, "V", 19, 1)
    cover = [line for line in lines if line.startswith("COVER")]
    assert summary["cover_without_total_support"] == 4 and all(line.endswith("ray 1 facet 1") for line in cover)
    assert any(line.startswith("COVER 0 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 0 0 0 ") for line in cover)
