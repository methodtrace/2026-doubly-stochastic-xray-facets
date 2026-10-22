# Facets of the X-ray polytope of doubly stochastic matrices: computational artifact

## Identifiers

- This version (1.0.1) of this artifact: https://doi.org/10.5281/zenodo.23440726
- All versions of this artifact: https://doi.org/10.5281/zenodo.23440725
- Preprint of the associated manuscript: arXiv:2610.30058, https://doi.org/10.48550/arXiv.2610.30058
- Published version: none yet
- Repository: https://github.com/methodtrace/2026-doubly-stochastic-xray-facets, tag v1.0.1
- Companion paper: arXiv:2610.24197

Code and data for every computational statement of the paper *Facets of the X-ray polytope of doubly stochastic
matrices* (authors and citation: `CITATION.cff`). The paper studies the polytope hat-P(k) of anti-diagonal X-rays of
k x k doubly stochastic matrices, the image of the Birkhoff polytope under the map that sums the entries along each
anti-diagonal (equivalently, the convex hull of the laws of i + sigma(i), sigma a permutation of {0, ..., k-1}, in
count coordinates). Its facets other than the coordinate facets come from potential inequalities I(u) with u a ray
of a tropical fan, the normal fan of a Minkowski sum of simplices. This artifact contains

- a certified facet test for the inequalities I(u), an exact and provably complete enumeration of the rays of the fan
  (the facet walk), and independent recomputations (a walk on the vertices of the Minkowski sum, in exact integer
  arithmetic and with cddlib; a line enumeration of the fan; exact facet lists from all permutation profiles for
  k <= 8);
- the complete ray lists for 5 <= k <= 12 and the certified classification of every ray (facet counts 2, 6, 13, 23,
  45, 85, 187, 426, 1033, 2633, 7041 for k = 2, ..., 12), with the census by step family, cores and heights;
- the checks of the facet criterion (Theorem 3.1 [thm:empty], Theorem 3.2 [thm:symdual], Theorem 3.3 [thm:F]) against
  the certified test, of the fiber of a facet (Theorem 3.9 [thm:fiber]) and the cubical fibers (Proposition 3.11
  [prop:cubical]), and of the facet ray without total support at k = 19 (Example 3.7 [ex:noTS], Remark 3.10
  [rem:samefacet]) with an exhaustive search over rays of height at most two;
- the checks of padding and cores (Theorem 4.2 [thm:pad]), grid refinement (Theorem 4.8 [thm:refine]) and its descent
  (Proposition 4.10 [prop:descent]), the non-step families (Proposition 4.11 [thm:fam]) and quadratic height
  (Theorem 4.13 [thm:height]);
- the checks behind the statements on step potentials (Theorem 5.13 [thm:S4], the three-step classification Theorem
  5.17 [thm:T3], nonincreasing tails Theorem 5.19 [thm:dec]) and on total support of staircases (Theorem 6.5 [thm:H]
  on Hall sets of block starts, Conjecture 6.2 [conj:T] on all staircases with k <= 24, statement E-flat on all tails
  with g_1 + ... + g_r <= 25, statement (S1) on all compositions with k <= 30);
- complete certificate streams with exact checkers for the computer-assisted theorems: statement E-flat for r <= 7
  steps (Theorem 6.8 [thm:Eflat6]) and statement (S1) for all real sequences with r <= 10 steps (Theorem 6.17
  [thm:S1cert]; r = 10 in a separate deposit file, see below), with the C generators that reproduce them.

Every computational claim of the paper, keyed by the number it carries in the current preprint and its LaTeX label
(for example "Theorem 3.3 [thm:F]", "Table SM1 [tab:cores]"; SM numbers are in the supplement), is listed in
`CLAIM_EVIDENCE_MAP.md` with the script that checks it, its data file and the tier that reaches the stated range.
The code itself cites results by label only.

## Requirements

Default (quick) tier, `python3 run.py`: CPython 3.12 or newer with `numpy`, `scipy` and `pytest`, and a C compiler
(`cc`, `gcc` or `clang`). The quick tier runs the test suite first and compiles `tools/c/attaining_structures.c`
(`enumerate_attaining_structures.py`) and `tools/c/search_total_support.c` (`search_total_support.py` and the
`--search` option of `check_facets_without_total_support.py`); `run.py` stops with a message naming what is missing.
The certificate tier (`--tier certificates`, below) needs only `numpy` and `scipy`: no pytest and no compiler.

Optional: `pycddlib` (package `pycddlib-standalone`, version 3 API) for the floating-point cross-checks with cddlib.
Without it the corresponding checks are skipped and reported as skipped (the paper tier passes `--require-cdd`, so
there a missing pycddlib is a failure). The other C programs (`boundary_classes.c`, `hplus_statement.c`,
`sign_cell_certificates.c`, `s1_certificate_tree.c`) are compiled only when a compiler is present or a medium, full or
paper tier option asks for them.

Tested versions: `pyproject.toml` pins numpy, scipy, pytest and pycddlib-standalone to the exact versions with which
every tier was last run (CPython 3.13, GCC 13.3). The quick tier has also passed with CPython 3.12, numpy 2.2 to 2.5
and scipy 1.16 to 1.18.

```sh
python3 -m pip install -e ".[test,cdd]"   # the pinned versions; or, unpinned:
python3 -m pip install numpy scipy pytest pycddlib-standalone
```

The scripts and tests run from a checkout without installation (they add `src/` to the import path). The checkers of
the certificate streams use only the standard library (the `.xz` streams are read with the `lzma` module).

## Audit path

Five commands, in this order, check the shipped files and the two main computational theorems (times measured on
2 cores of a 2.1 GHz class x86_64 server, CPython 3.13):

```sh
sha256sum -c --quiet MANIFEST.sha256                         # 1. every shipped file is intact (prints nothing; seconds)
python3 scripts/reproduce_facet_counts.py                     # 2. Table 2 [tab:counts]: certified classification of all deposited rays k <= 12, facet walk k <= 9, profiles k <= 8 (about 1 minute)
python3 scripts/check_cells_depth7.py --full                  # 3. Theorem 6.8 [thm:Eflat6]: E-flat sign-cell certificate streams r <= 7 (about 3 minutes)
python3 scripts/check_s1_certificates.py --full               # 4. Theorem 6.17 [thm:S1cert], r <= 9: (S1) certificate trees (about 7 minutes)
python3 scripts/check_s1_depth10.py --deposit s1_depth10_certificates.zip --workers 2   # 5. Theorem 6.17 [thm:S1cert], r = 10 (about 35 minutes; needs the separate deposit file)
```

Without the deposit file, command 5 checks only the r = 10 index (`python3 scripts/check_s1_depth10.py`, seconds);
the r = 10 claim is then "index only" (see the status labels in `CLAIM_EVIDENCE_MAP.md`). On macOS use
`shasum -a 256 -c MANIFEST.sha256` for command 1. Commands 1 to 5 (command 2 without the facet walk and the profiles,
commands 3 and 4 without rerunning the C generators) are what `python3 run.py --tier certificates` runs.

## Running

```sh
python3 run.py                         # quick tier: tests + every script at its default range
python3 run.py --tier certificates     # verify the shipped certificates only, recompute nothing (see below)
python3 run.py --tier medium           # larger ranges, a few hours on 2 cores
python3 run.py --tier full             # every script with --full (about 150 core-hours)
python3 run.py --tier paper            # every range stated in the paper, as independent pieces (next section)
python3 run.py --only cells,padding --tier full     # selected steps
```

Quick tier runtime: under 25 minutes on 2 cores; the test suite alone takes about 3 to 7 minutes.

The certificate tier verifies what is shipped without recomputing it: the SHA-256 of every file in
`MANIFEST.sha256`; the certified facet test on every deposited ray k <= 12 against the deposited outcomes and
certificates (no facet walk beyond k = 4, no line enumeration); every record of the E-flat sign-cell streams r <= 7
and of the (S1) certificate trees r <= 9 (the C generators are not rerun); and the r = 10 domain index, and every
r = 10 domain tree when `--deposit s1_depth10_certificates.zip` is given. It needs no C compiler and no pytest. On
the machine above the (S1) trees r <= 9 take 408 s and the E-flat streams r <= 7 take 183 s, about 10 minutes in all.
Steps that check less than the paper's range are marked in the runner summary (for example "index only" for r = 10
without the deposit file).

Every step writes `results/<step>/<name>.json` with its counts and a pass/fail entry per check, and `run.py` writes
`results/run_summary.json`. Each script can also be run alone; `python3 scripts/<script>.py --help` lists its
options, every script takes `--out DIR`, small default ranges, and `--full` for its full-tier range.

## The paper tier

`python3 run.py --tier paper` runs every range stated in the paper, including the facet walk at k = 12, the line
enumeration at k = 10, the vertex walks at k = 11 and 12, the staircase sweep up to k = 24 and the regeneration of
every certificate tree (the E-flat stream r = 7, the (S1) trees r <= 9 and the 441 domain trees r = 10, from which
it builds the deposit file `s1_depth10_certificates.zip`). The tier is a list of independent pieces:

```sh
python3 run.py --tier paper --list               # the pieces as JSON: id, step, arguments, dependencies
python3 run.py --tier paper --piece ID[,ID...]   # run selected pieces
python3 run.py --tier paper --workers 8          # run every piece that has not passed yet
python3 run.py --tier paper --out DIR            # results under DIR instead of results/
```

Each piece writes its results, its log (`log.txt`) and `piece.json` to `OUT/pieces/<id>/`; `piece.json` records the
status, the wall time, the CPU time and the peak memory of the largest process, with the number of workers and a
description of the machine. A rerun skips every piece whose `piece.json` reports a pass, so an interrupted run resumes
where it stopped. Pieces can run on different machines and in any order that respects the dependencies listed by
`--list` (the merge of the staircase pieces needs every staircase piece, and the deposit file needs every r = 10
regeneration piece); copy the finished `OUT/pieces/<id>/` folders into one `OUT/pieces/` before running a piece
that depends on them. The result of a piece does not depend on the machine or on the number of workers.

Measured cost of the paper tier:

| Machine | Cores | Wall time | CPU time |
|---|---|---|---|
| Linux x86_64 cluster nodes, four lanes run at the same time | 4 x 16 | about 6 hours elapsed | 80 hours |
| Apple silicon desktop, arm64 (4 performance, 6 efficiency, 2 super cores), 16 GB | 12 (6 for the four heavy pieces) | about 8 hours of running time, run in budgeted sessions | 47 hours |

Largest process: 4.6 GB on the cluster, 3.7 GB on the Mac. The two runs agree file for file (uncompressed SHA-256) and
number for number. The heaviest pieces are the staircase sweep at k = 23, 24 and the walks at k = 12.

## Reproducibility standard

Two runs of a piece, on any two machines, must agree as follows:

- every file the piece regenerates (ray lists, facet lists, certificate streams and trees, the deposit file) has the
  same SHA-256 of its uncompressed content; compressed bytes depend on the version of the zlib or xz library and are
  not compared;
- every number in its result files is identical;
- the floating-point cross-checks (the cddlib vertex walk and facet lists, and the floating-point LP comparison of
  the elimination routine of the cell enumeration) only have to pass.

`tools/fingerprint.py OUT` records these fingerprints for a paper-tier run, and `tools/fingerprint.py --compare
A.json B.json` compares two runs.

## Scripts

| Script | Reproduces | Quick tier | Full tier (2.1 GHz class core) |
|---|---|---|---|
| `reproduce_facet_counts.py` | Table 2 [tab:counts] for k = 2..12: certified classification of the deposited rays (k = 5..12), exact facet lists from all permutation profiles (k = 2..8), facet iff total support on all rays, the facet walk, line enumeration | classification k <= 12, profiles k <= 8, walk k <= 9 (2 min) | walk k <= 11 (45 min), line enumeration k <= 9; `--walk-kmin 12 --walk-kmax 12`: 9 core-hours; `--line-kmin 10 --line-kmax 10`: 8 core-hours |
| `reproduce_census.py` | Table 3 [tab:census] | k <= 12 (seconds) | same |
| `reproduce_shapes.py` | Table SM2 [tab:shapes], hinges (Proposition 5.8 [prop:hinges]), heights (Remark 2.5 [rem:multiplicity]) | k <= 12 (seconds) | same |
| `crosscheck_vertex_walk.py` | Table 4 [tab:vw]: independent vertex walk (exact and cddlib), cddlib facet lists with exact rows | k <= 9 (1 min) | exact k <= 11 (1 hour); `--cdd-kmin 11 --cdd-kmax 12`: 6 core-hours; `--kmin 12 --kmax 12`: exact walk at k = 12 (about 7.5 core-hours; 29 minutes on 16 cores in the reproduction run) |
| `check_staircase_conjectures.py` | Conjecture 6.2 [conj:T] on all staircases (Table SM3 [tab:evidence]), the rigid-corner certificate, Theorem 5.13 [thm:S4] counts, Lemma 5.4 [lem:last], Lemma 5.5 [lem:starts]; runs in pieces (`--k`, `--r`, `--merge`) | k <= 13 (1 min) | k <= 24: 100 to 150 core-hours |
| `check_nonincreasing_tails.py` | Theorem 5.19 [thm:dec], Corollary 5.20 [cor:dec] | k <= 14 (10 s) | k <= 20 (2 min on 2 cores) |
| `check_three_step_theorem.py` | Theorem 5.17 [thm:T3]: proof steps and statement | k <= 14 / 22 (1 min) | k <= 30 / 40: about 3 core-hours |
| `check_hall_theorem.py` | Theorem 6.5 [thm:H], its proof inequalities | k <= 13 (seconds) | k <= 18: about 30 minutes |
| `check_boundary_classes.py` | Corollary 6.16 [cor:probB] ((S1) totals), Lemma 6.14 [lem:R] | Python k <= 16 (20 s) | C program k <= 30: about 40 minutes |
| `check_cells.py` | Theorem 6.8 [thm:Eflat6] for r <= 6, Section SM6.7 [app:cells] (Table SM4 [tab:cells], checks (i), (ii)), E-flat on all tails with sum <= 25 | r <= 5 (30 s) | r <= 6 and sums <= 25: about 1 hour |
| `check_cells_depth7.py` | Theorem 6.8 [thm:Eflat6] (r <= 7), Corollary 6.11 [cor:probA], Lemma 6.10 [lem:cellpair], Table SM4 [tab:cells]: sign-cell certificate streams, a second (matching-based) verdict on every cell, and the regeneration of the streams by the C generator (r = 7 with `--full`) | r <= 6 (10 s) | r = 7 (4 minutes; regeneration of the r = 7 stream with a C compiler) |
| `check_staircase_subclasses.py` | Proposition 6.12 [prop:Eflatclasses], Proposition 6.13 [prop:stairbip]: unimodal tails, unique maximizing pairs, superadditive prefixes, bipartite allowed component; a further tail not cited in the paper | r <= 7, k <= 10 (3 s) | r <= 9, k <= 12 (20 s) |
| `enumerate_attaining_structures.py` | (S1) for every sequence with r <= 10 by an exact enumeration of attaining structures (a second route to Theorem 6.17 [thm:S1cert]); catalogue of attaining structures r <= 6 | r <= 7 (5 s) | r <= 9 (8 min); `--r10-chunks 6`: r = 10, 2.3 core-hours |
| `check_s1_certificates.py` | Theorem 6.17 [thm:S1cert] (r <= 9): (S1) certificate trees; with a C compiler the generator regenerates every tree r <= rmax and one r = 10 domain byte for byte | r <= 7 (5 s) | r <= 9 (7 minutes) |
| `check_s1_depth10.py` | Theorem 6.17 [thm:S1cert] (r = 10): domain index; with `--deposit`, the 441 domain trees of the deposit file (files identified by the SHA-256 of their uncompressed tree); `--make-deposit` builds the deposit file | index (seconds) | same; `--deposit`: about 1.5 core-hours |
| `regenerate_s1_depth10.py` | Theorem 6.17 [thm:S1cert] (r = 10): regenerates domain trees with the C generator, checks the SHA-256 of each uncompressed tree against the index and verifies each tree completely; `--domains A:B` selects a range of domains | the 20 smallest domains (seconds) | all 441 domains: a few core-hours (paper tier pieces `s1-regenerate-NN`) |
| `check_facet_criterion.py` | Theorem 3.1 [thm:empty], Theorem 3.2 [thm:symdual], Theorem 3.3 [thm:F], Proposition 3.4 [prop:raybip], Corollary 3.6 [cor:facetbip], on random potentials (integer potentials and exact integer Hungarian duals) and on all rays | random k <= 14, rays k <= 9 (15 s) | random k <= 16 (minutes) |
| `check_facets_without_total_support.py` | Example 3.7 [ex:noTS], Remark 3.10 [rem:samefacet] (with `--search`: rays of height <= 2) | search k <= 12 (5 s) | search k <= 18 (28 core-minutes) |
| `search_total_support.py` | exhaustive searches for failures of total support (Remark 3.10 [rem:samefacet] and data) | 15 s | values (13, 2), (16, 2), increments k <= 16 (minutes) |
| `check_facet_fibers.py` | Theorem 3.9 [thm:fiber], Proposition 3.11 [prop:cubical]: facet fibers, cubical fibers on 19f points (the script calls f m) | f <= 2 (5 s) | f <= 3 (30 s); `--mmax 6`: f <= 6 |
| `check_height_two_rays.py` | Remark 3.12 [rem:padcount]: lower bound for the number of height-two rays | k <= 24 (2 s) | k <= 40 (30 s) |
| `check_padding.py` | Theorem 4.2 [thm:pad], Corollary 4.3 [cor:pad], Corollary 4.4 [cor:cores], Table SM1 [tab:cores], block lifts of the cores (not cited in the paper) | k <= 12 (15 s) | adds proof steps k = 11, paddings to k = 13, block lifts k = 12 |
| `check_refinement.py` | Theorem 4.8 [thm:refine], Corollary 4.9 [cor:finer], Lemma 4.5 [lem:sprtight], Lemma 4.6 [lem:blockgraph], Lemma 4.7 [lem:flatlevel] | facet rays k <= 8, m = 2, 3 (25 s) | k <= 9 (1.5 min) |
| `check_refinement_descent.py` | Proposition 4.10 [prop:descent]: descent under refinement and its counterexample; ray-only padding (Theorem 4.2(b) [thm:pad]) | k <= 5 (2 s) | k <= 6, padding k <= 7 (5 s) |
| `check_families.py` | Proposition 4.11(a)-(c), (e) [thm:fam], Section SM4.5 [sec:SM-famproof] and the coverage after Proposition 4.11 | k <= 16 / 24 (seconds) | k <= 30 (minutes) |
| `check_family_distinctness.py` | Proposition 4.11(d) [thm:fam], the patterns of Lemma 3.8 [lem:distinct] | k <= 24 (seconds) | k <= 40 (10 s) |
| `check_unbounded_height.py` | Theorem 4.13 [thm:height] for a = 1 (height b on 2b + 2 points) and the bound floor((k - 2)/2); the certified facet test on G_ab | b <= 10 (seconds) | b <= 40 (2 min) |
| `check_quadratic_family.py` | Theorem 4.13 [thm:height]: quadratic family G_ab (facet ray iff a <= b + 1), height bound floor((k - 2)/3)^2 | b <= 8, k <= 20 (5 s) | b <= 12, k <= 30 (10 s) |
| `check_two_point_lifts.py` | Problem 6.1 [prob:C] and Section SM6.1 [sec:s_census-lifts] (two-point lifts, a rule verified by exact computation): 962 of the 2,808 cores at k = 12 | grid on 3 parents (20 s) | grid on 8 parents (20 s) |
| `check_binary_increments.py` | binary interval-cancellation certificate for rays (not cited in the paper) | k <= 10 (1 s) | k <= 12 (3 s) |
| `check_examples.py` | worked examples of Sections 2 to 6 and of the supplement, Remark 2.5 [rem:multiplicity] (the tail-sum reading with repeated thresholds), Remark 4.12 [prop:P] | seconds | same |

Times are indicative. The facet walk at k = 12 keeps its state in its output folder and resumes after an
interruption; within the paper tier, every piece resumes at the level of pieces. All other steps simply rerun.

## Exactness

Every count and every verdict is computed in exact integer or rational arithmetic, with certificates:

- the facet test reports for each potential one of: no tight permutation; two anti-diagonals without allowed edge
  (one suffices by Theorem 3.1(b) [thm:empty]; the test uses two for every potential, staircases included); an
  equation; a rank modulo the prime 2^31 - 1 reaching the facet dimension (a modular rank never exceeds the rational
  rank, and the incidence part is totally unimodular); or an exact rational recomputation;
- the facet walk computes ranks modulo 2^30 - 35, exact for k <= 25 by Hadamard's inequality, and checks that every
  ridge is found from both sides; its completeness is Proposition A.1 [prop:complete];
- the cell enumeration decides the nonemptiness of each cell by exact Fourier-Motzkin elimination;
- the facet criterion uses exact integer ranks, and its random test potentials include exact integer optimal duals
  of assignment problems (Hungarian method); the enumeration of attaining structures backs every feasibility
  decision by an integer witness and every infeasibility by a Gordan multiplier, both checked exactly in the C
  program, and the deposited catalogue is rechecked by independent Python code;
- the C search over potentials (`tools/c/search_total_support.c`) computes ranks modulo two primes, as stated in
  Remark 3.10 [rem:samefacet]: a ray or facet it reports is certified (and rechecked exactly in Python), and a ray or
  facet could be missed only if both modular ranks fell below the rational rank;
- the certificate streams (E-flat for r <= 7, (S1) for r <= 10) are checked record by record in unbounded integer
  arithmetic: every rejected branch carries a Gordan multiplier or a strict row in the span of the equalities, every
  realized branch an integer witness or a combinatorial cut, and the terminal intervals must partition the whole
  tree; the checkers share no code with the C generators that produced the streams
  (`tools/c/sign_cell_certificates.c`, `tools/c/s1_certificate_tree.c`, which reproduce the deposited streams byte
  for byte); `.xz` streams are also checked against the SHA-256 of their uncompressed content.

Floating point enters only in cross-checks that are labeled as such: the cddlib vertex walk and cddlib facet lists
(every reported ray and row is verified exactly, and the wall detection of the cddlib vertex walk uses a
floating-point rank), and the cross-check of the elimination routine of the cell enumeration against a linear
program (scipy).

## Contents

- `src/dsxray/`: the package (one module per topic; see the module docstrings and `src/dsxray/__init__.py`).
- `scripts/`: one script per group of claims, named by what it reproduces; `scripts/common.py` holds shared helpers.
  The docstring of each script lists the claims it checks by their LaTeX labels.
- `tests/`: unit and regression tests (`python3 -m pytest -q`).
- `run.py`: the tiers and the pieces of the paper tier.
- `tools/c/`: `boundary_classes.c` (statement (S1), exhaustive over compositions), `attaining_structures.c`
  ((S1) for all sequences with r <= 10 steps, and the catalogue of attaining structures), `search_total_support.c`
  (exhaustive searches over potentials with bounded values or height), `hplus_statement.c` (statement (H+), not used
  in the paper), and the generators of the certificate streams `sign_cell_certificates.c` (E-flat, Theorem 6.8
  [thm:Eflat6]) and `s1_certificate_tree.c` ((S1), Theorem 6.17 [thm:S1cert]), compiled on demand by the scripts.
- `tools/fingerprint.py`: fingerprint of a paper-tier run and comparison of two runs (see "Reproducibility
  standard").
- `data/rays/rays_kNN.json`: every ray of the fan for k = 5..12 (normalized u(0) = u(k-1) = 0, primitive integer).
- `data/facets/facets_kNN.jsonl` (k = 8..11) and `facets_k12.jsonl.gz`: the certified classification of every ray
  (one JSON object per ray: ray, dim Z, facet, certificate, total support, edges outside every tight permutation,
  empty anti-diagonals).
- `data/cells/cells_rN.json` (r = 1..5) and `cells_r6.json.gz`: the nonempty cells of the E-flat enumeration with
  their verdicts.
- `data/cells/sign_certificates_rN.txt.gz` (r = 1..6) and `sign_certificates_r7.txt.xz` (38 MB; 413 MB
  uncompressed): complete sign-cell certificate streams for statement E-flat (format in
  `src/dsxray/sign_certificates.py`); `sign_certificates_r7.json`: counts and the SHA-256 of the uncompressed stream.
- `data/cores/cores_kNN.json`: the cores of the facet rays (increment patterns) for k = 3..12;
  `data/cores/block_lifts_kNN.json` (k = 5..12): each core with its block lifts (the block faces of hat-P(k),
  isomorphic to hat-P(k - b), that its facet meets in a facet; data only, not cited in the paper).
- `data/counterexamples/potentials.json`: the potentials of Example 3.7 [ex:noTS] and Remark 3.10 [rem:samefacet]
  (k = 19), two potentials at k = 17 and the family R(m) of facet rays without total support (k = 22..41).
- `data/attaining_structures/structures_rN.jsonl` (r = 1..5) and `structures_r6.jsonl.gz`: every realizable
  attaining structure of a sequence E_0..E_r (pairs attaining at each level s = 0..2r-1) with an integer witness
  (E_0 = E_1 = 0).
- `data/families/quadratic_height_family.json`: the certified facet test on the quadratic family G_ab for b <= 8.
- `data/s1/s1_certificates_rN.txt.gz` (r = 1..8) and `s1_certificates_r9.txt.xz` (34 MB; 400 MB uncompressed):
  certificate trees for statement (S1) (format in `src/dsxray/s1_certificates.py`); `s1_certificates.json`: terminal
  counts and SHA-256 of each uncompressed tree. `s1_depth10_domains.json`, `s1_depth10_certificates.sha256` and
  `s1_depth10_trees.sha256`: index and checksums of the 441 domain trees for r = 10, which are in the separate deposit
  file (next section).
- `data/summary/`: counts from the full-tier runs (facet counts and census, vertex walk table, staircase sweeps,
  Theorem 6.5 [thm:H], (S1) and (H+) exhaustive counts, E-flat). Scripts compare their output with these files where
  the ranges overlap.
- `MANIFEST.sha256`: SHA-256 checksums of all files (not of the separate deposit file).

## The r = 10 certificates for (S1): separate deposit file

The certificate tree for (S1) with r = 10 steps is too large for this archive (441 domain files, 355 MB compressed,
4.7 GB uncompressed). It is deposited as the separate file `s1_depth10_certificates.zip` (in the Zenodo record and attached to the
GitHub release v1.0.1), which
contains the directory `s1_depth10_certificates/` with one file `<prefix>.txt.xz` per domain (the prefix is the mask
prefix of length 6, for example `1_1_3_2_7_5.txt.xz`), together with copies of `data/s1/s1_depth10_domains.json`,
`s1_depth10_certificates.sha256` and `s1_depth10_trees.sha256`. This archive holds the index, the checksums, the
generator and the checker. To verify:

```sh
python3 scripts/check_s1_depth10.py                                              # the index (seconds)
python3 scripts/check_s1_depth10.py --deposit s1_depth10_certificates.zip --domains 320 # hashes + 320 smallest domains (1 minute)
python3 scripts/check_s1_depth10.py --deposit s1_depth10_certificates.zip --workers 4   # every domain (about 1.5 core-hours)
```

`--deposit` also accepts the unpacked directory `s1_depth10_certificates/`; there, `sha256sum -c
s1_depth10_certificates.sha256` checks the files. The full verification checks every terminal record of every
domain exactly (as `check_s1_certificates.py` does for r <= 9), that each domain tree covers exactly its interval,
the SHA-256 of each uncompressed tree, and the counts of the index; the index check confirms that the 441 intervals
partition the whole tree [0, 572834992892439375). The deposit file can also be rebuilt from scratch:
`regenerate_s1_depth10.py` regenerates the domain trees with the C generator (in pieces with `--domains A:B`), and
`check_s1_depth10.py --make-deposit DIR OUT.zip` packs them; the paper tier does both. A rebuilt file matches the
deposited one in the SHA-256 of every uncompressed tree; its compressed bytes may differ (they depend on the version of
the xz library).

## License

Code (`src/`, `scripts/`, `tests/`, `tools/`, `run.py`): Apache License 2.0 (`LICENSE`).
Data (`data/`): Creative Commons Attribution 4.0 International (`LICENSE-DATA`).
