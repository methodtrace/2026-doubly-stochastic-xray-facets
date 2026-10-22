"""Attaining structures: C enumeration against the deposited catalogue, independent witness check, (S1)."""
import os
import shutil
import subprocess

import pytest

from common import DATA, ROOT, read_jsonl
from dsxray import attaining_structures as AS


def catalogue(r):
    p = os.path.join(DATA, "attaining_structures", f"structures_r{r}.jsonl" + (".gz" if r >= 6 else ""))
    return [(rec["levels"], rec["witness"]) for rec in read_jsonl(p)]


@pytest.mark.parametrize("r,n", [(1, 1), (2, 3), (3, 13), (4, 75), (5, 555), (6, 4421)])
def test_deposited_catalogue(r, n):
    chk = AS.check_catalogue(catalogue(r), r)
    assert chk["records"] == chk["distinct"] == n and chk["witness_failures"] == 0 and chk["s1_failures"] == 0


def test_brute_force_subset_r4():
    cat = {tuple(tuple(tuple(p) for p in lev) for lev in s) for s, _ in catalogue(4)}
    assert AS.brute_force_structures(4, 5) <= cat


@pytest.fixture(scope="module")
def exe(tmp_path_factory):
    cc = shutil.which("cc") or shutil.which("gcc")
    if not cc:
        pytest.skip("no C compiler")
    out = str(tmp_path_factory.mktemp("bin") / "attaining_structures")
    subprocess.run([cc, "-O2", "-o", out, os.path.join(ROOT, "tools", "c", "attaining_structures.c")], check=True)
    return out


def test_c_catalogue_equals_deposited(exe):
    for r in (3, 4, 5):
        lines = subprocess.run([exe, "C", str(r)], capture_output=True, text=True, check=True).stdout.splitlines()
        got = {AS.parse_c_line(x)[0] for x in lines if x.startswith("AS ")}
        assert got == {tuple(tuple(tuple(p) for p in lev) for lev in s) for s, _ in catalogue(r)}


def test_c_s1_small(exe):
    import json
    for r in range(1, 7):
        out = subprocess.run([exe, "S", str(r)], capture_output=True, text=True, check=True).stdout.splitlines()
        assert json.loads(out[-1])["s1_failures"] == 0 and len(out) == 1
