/* Statement (S1) for staircases, exhaustive, in exact integer arithmetic.
   (S1): if the boundary steps 0..r-1 form at least two classes, some q in 0..r-1 has at most one attaining
   partner l in 0..r, i.e. #{l : E_q + E_l = M(q+l)}  <= 1, where M(s) = max_{p+q=s} E_p + E_q.
   Also counts the cases where such a q lies outside the class of 0.
   Mode A: all compositions of k with 1 <= r <= k-2, for KMIN <= k <= KMAX.         usage: boundary_classes A KMIN KMAX
   Mode P: pattern lifts E_p = P_p + (V+1) p, P in {0..V}^(r+1), min P = 0, 2 <= r <= RMAX
           (by affine invariance these are end sequences of compositions).        usage: boundary_classes P V RMAX
   Build: cc -O2 -o boundary_classes boundary_classes.c */
#include <stdio.h>
#include <stdlib.h>

static int par[128];
static int find(int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; }
static long long nmulti, nfail, noutside;

static void test(const long long *E, int r) {
  long long M[256]; static int att[64][64];
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
  int c0 = find(0), multi = 0;
  for (int p = 0; p < r; p++) if (find(p) != c0) multi = 1;
  if (!multi) return;
  nmulti++;
  int found = 0, outside = 0;
  for (int q = 0; q < r; q++) {
    int c = 0;
    for (int l = 0; l <= r; l++) c += att[q][l];
    if (c <= 1) { found = 1; if (find(q) != c0) outside = 1; }
  }
  if (outside) noutside++;
  if (!found) {
    nfail++;
    if (nfail <= 10) { printf("  failure E="); for (int i = 0; i <= r; i++) printf("%lld,", E[i]); printf("\n"); }
  }
}

int main(int argc, char **argv) {
  if (argc != 4) { fprintf(stderr, "usage: boundary_classes A KMIN KMAX | boundary_classes P V RMAX\n"); return 2; }
  if (argv[1][0] == 'A') {
    int kmin = atoi(argv[2]), kmax = atoi(argv[3]);
    for (int k = kmin; k <= kmax; k++) {
      nmulti = nfail = noutside = 0; long long n = 0;
      for (long long cuts = 0; cuts < (1LL << (k - 1)); cuts++) {
        long long E[64]; int r = 0;
        for (int i = 0; i < k - 1; i++) if (cuts >> i & 1) E[r++] = i;
        E[r] = k - 1;
        if (r < 1 || r > k - 2) continue;
        n++; test(E, r);
      }
      printf("{\"k\": %d, \"compositions\": %lld, \"several_classes\": %lld, \"solitary_outside_class_of_0\": %lld, \"s1_failures\": %lld}\n",
             k, n, nmulti, noutside, nfail);
      fflush(stdout);
    }
  } else {
    int V = atoi(argv[2]), RM = atoi(argv[3]);
    for (int r = 2; r <= RM; r++) {
      nmulti = nfail = noutside = 0;
      long long tot = 1;
      for (int i = 0; i <= r; i++) tot *= (V + 1);
      for (long long code = 0; code < tot; code++) {
        int P[64]; long long c = code;
        for (int i = 0; i <= r; i++) { P[i] = c % (V + 1); c /= (V + 1); }
        int mn = P[0];
        for (int i = 1; i <= r; i++) if (P[i] < mn) mn = P[i];
        if (mn != 0) continue;
        long long E[64];
        for (int p = 0; p <= r; p++) E[p] = P[p] + (long long)(V + 1) * p;
        test(E, r);
      }
      printf("{\"V\": %d, \"r\": %d, \"several_classes\": %lld, \"s1_failures\": %lld}\n", V, r, nmulti, nfail);
      fflush(stdout);
    }
  }
  return 0;
}
