/* Statement (H+) on all compositions of k with 1 <= r <= k-2, KMIN <= k <= KMAX (exact integer arithmetic).
   Verdicts: connected (one class of boundary steps), deficient (a class K with |A(K)| < |K|), tight_set (a balanced
   class with a tight edge from an end E_l, l in A(K), to an index outside the starts of K), FAIL (neither).
   Tightness by the corner lemma: E_l is tight with index j of block p iff max(S_p, M(p+l-1)+1-E_l) <= j <= E_p
   (for (p,l) = (0,0): 0 <= j <= E_0); j lies among the starts of K iff j = S_p with p-1 in K.
   usage: hplus_statement KMIN KMAX        build: cc -O2 -o hplus_statement hplus_statement.c */
#include <stdio.h>
#include <stdlib.h>

static int par[128];
static int find(int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; }

static int verdict(const long long *E, int r) {
  long long M[256]; int att[64][64];
  for (int s = 0; s <= 2 * r; s++) {
    long long m = -(1LL << 60);
    for (int p = (s - r > 0 ? s - r : 0); p <= (s < r ? s : r); p++) { long long v = E[p] + E[s - p]; if (v > m) m = v; }
    M[s] = m;
  }
  for (int p = 0; p <= r; p++) for (int q = 0; q <= r; q++) att[p][q] = (E[p] + E[q] == M[p + q]);
  for (int p = 0; p < r; p++) par[p] = p;
  for (int s = 0; s <= 2 * r; s++) {
    int first = -1;
    for (int p = (s - r > 0 ? s - r : 0); p <= (s < r ? s : r); p++)
      if (p < r && att[p][s - p]) { if (first < 0) first = p; else par[find(p)] = find(first); }
  }
  int root[64], nc = 0;
  for (int p = 0; p < r; p++) if (find(p) == p) root[nc++] = p;
  if (nc == 1) return 0;
  long long S[64]; S[0] = 0;
  for (int p = 1; p <= r; p++) S[p] = E[p - 1] + 1;
  int anyleak = 0;
  for (int c = 0; c < nc; c++) {
    int inK[64] = {0}, nK = 0, inA[64] = {0}, nA = 0;
    for (int q = 0; q < r; q++) if (find(q) == root[c]) { inK[q] = 1; nK++; }
    for (int q = 0; q < r; q++) if (inK[q]) for (int l = 0; l <= r; l++) if (att[q][l] && !inA[l]) { inA[l] = 1; nA++; }
    if (nA < nK) return 1;
    if (nA == nK) {
      for (int l = 0; l <= r && !anyleak; l++) if (inA[l]) for (int p = 0; p <= r; p++) {
        long long lo = (p == 0 && l == 0) ? 0 : S[p];
        if (p + l > 0) { long long t = M[p + l - 1] + 1 - E[l]; if (t > lo) lo = t; }
        if (lo > E[p]) continue;
        if (lo == E[p] && lo == S[p] && p >= 1 && inK[p - 1]) continue;
        anyleak = 1; break;
      }
    }
  }
  return anyleak ? 2 : 3;
}

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: hplus_statement KMIN KMAX\n"); return 2; }
  int kmin = atoi(argv[1]), kmax = atoi(argv[2]);
  for (int k = kmin; k <= kmax; k++) {
    long long cnt[4] = {0}, n = 0;
    for (long long cuts = 0; cuts < (1LL << (k - 1)); cuts++) {
      long long E[64]; int r = 0;
      for (int i = 0; i < k - 1; i++) if (cuts >> i & 1) E[r++] = i;
      E[r] = k - 1;
      if (r < 1 || r > k - 2) continue;
      n++; int v = verdict(E, r); cnt[v]++;
      if (v == 3 && cnt[3] <= 10) { printf("  failure k=%d E=", k); for (int i = 0; i <= r; i++) printf("%lld,", E[i]); printf("\n"); }
    }
    printf("{\"k\": %d, \"compositions\": %lld, \"connected\": %lld, \"deficient\": %lld, \"tight_set\": %lld, \"FAIL\": %lld}\n",
           k, n, cnt[0], cnt[1], cnt[2], cnt[3]);
    fflush(stdout);
  }
  return 0;
}
