"""Certificate streams for E-flat (sign cells) and (S1) (attaining-set trees): acceptance and rejection."""
import gzip
import os

import pytest

from dsxray import attaining_structures as AS
from dsxray import cell_enumeration as CE
from dsxray import s1_certificates as S1
from dsxray import sign_certificates as SC

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def read_lines(path):
    with gzip.open(path, "rt") as f:
        return f.read().splitlines()


def write_lines(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def test_forms_and_verdicts():
    assert [len(SC.comparison_forms(r)) for r in range(1, 8)] == [0, 1, 3, 7, 13, 22, 34]
    for t in [(3, 1, 2, 2), (1, 1), (2, 5, 3, 1, 1), (4, 5, 4, 1, 5, 5), (2, 2, 3, 1, 1, 4, 2)]:
        assert SC.verdict_of_tail(t) == CE.eflat_tail(t)


def test_sign_stream_r5_matches_deposited_cells():
    from common import read_json
    dep = read_json(os.path.join(DATA, "cells", "cells_r5.json"))
    s = SC.verify_stream(os.path.join(DATA, "cells", "sign_certificates_r5.txt.gz"), 5, dep, CE.eflat_tail)
    assert s["cells"] == 2298 and s["cell_verdicts"] == {"nohall": 586, "ts": 434, "ok": 1278}


@pytest.mark.parametrize("kind", ["gordan", "witness", "drop", "summary"])
def test_sign_stream_corruptions_rejected(tmp_path, kind):
    lines = read_lines(os.path.join(DATA, "cells", "sign_certificates_r4.txt.gz"))
    if kind == "gordan":
        i = next(i for i, x in enumerate(lines) if " G " in x)
        f = lines[i].split(); f[4] = f[4].split(":")[0] + ":1/7"; lines[i] = " ".join(f)
    elif kind == "witness":
        i = next(i for i, x in enumerate(lines) if x.startswith("L ") and "0" not in x.split()[2])
        f = lines[i].split(); f[4], f[5] = f[5], f[4]; lines[i] = " ".join(f)
    elif kind == "drop":
        del lines[len(lines) // 2]
    else:
        lines = lines[:-1]
    p = tmp_path / "s.txt"; write_lines(p, lines)
    with pytest.raises((AssertionError, ValueError)):
        SC.verify_stream(str(p), 4)


def test_s1_tree_r6_with_probes():
    probes = [[0, 0, 1, 0, 2, 1, 0], [0, 0, -1, 1, 0, 0, 2], [0, 0, 0, 0, 0, 0, 0], [0, 0, 3, 1, 1, 3, 0]]
    s = S1.verify_tree(os.path.join(DATA, "s1", "s1_certificates_r6.txt.gz"), 6, (), probes)
    assert s["terminal_records"] == 9293 and s["probes_checked"] == 4
    assert all(AS.s1_holds(AS.structure_of(E), 6) for E in probes)


@pytest.mark.parametrize("kind", ["gordan", "cut", "drop", "summary"])
def test_s1_tree_corruptions_rejected(tmp_path, kind):
    lines = read_lines(os.path.join(DATA, "s1", "s1_certificates_r5.txt.gz"))
    if kind == "gordan":
        i = next(i for i, x in enumerate(lines) if x.startswith("G "))
        f = lines[i].split(); f[2] = f[2].split(":")[0] + ":1/9"; lines[i] = " ".join(f)
    elif kind == "cut":
        i = next(i for i, x in enumerate(lines) if x.startswith("Z "))
        lines[i] = "C " + lines[i].split()[1]
    elif kind == "drop":
        del lines[len(lines) // 2]
    else:
        lines = lines[:-1]
    p = tmp_path / "t.txt"; write_lines(p, lines)
    with pytest.raises(AssertionError):
        S1.verify_tree(str(p), 5)


def test_depth10_index_shape():
    from common import read_json
    idx = read_json(os.path.join(DATA, "s1", "s1_depth10_domains.json"))
    groups, radices, widths = S1.tree_shape(10)
    assert radices[:6] == [1, 1, 3, 3, 7, 7] and widths[0] == idx["tree_interval"][1]
    assert len(idx["domain_list"]) == 441


def test_generators_reproduce_streams(tmp_path):
    import hashlib
    import shutil
    import subprocess
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if not cc:
        pytest.skip("no C compiler")
    tools = os.path.join(ROOT, "tools", "c")
    s1 = str(tmp_path / "s1"); cells = str(tmp_path / "cells")
    subprocess.run([cc, "-O2", "-o", s1, os.path.join(tools, "s1_certificate_tree.c")], check=True)
    subprocess.run([cc, "-O2", "-o", cells, os.path.join(tools, "sign_cell_certificates.c")], check=True)
    from common import read_json
    meta = read_json(os.path.join(DATA, "s1", "s1_certificates.json"))
    out = subprocess.run([s1, "S", "5", "-1", "0", "1", "1", "1", "0", "_"], capture_output=True, check=True).stdout
    assert hashlib.sha256(out).hexdigest() == meta["5"]["uncompressed_sha256"]
    out = subprocess.run([cells, "4", "1", "0", "1"], capture_output=True, check=True).stdout
    with gzip.open(os.path.join(DATA, "cells", "sign_certificates_r4.txt.gz"), "rb") as f:
        assert out == f.read()
