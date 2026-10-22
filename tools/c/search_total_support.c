/* Exhaustive search over integer potentials for failures of total support.

   Usage: search_total_support V k H          (all potentials with values in {0..H}, u(0) = 0)
          search_total_support I k H [LAST]   (all normalized potentials of height <= H)
   Build: cc -O2 -o search_total_support search_total_support.c

   For a potential u on {0..k-1}: v(t) = min_{i+j=t} u(i)+u(j), tight graph E (pairs attaining v), a tight
   permutation by augmenting paths, and the allowed graph A (tight edges in some tight permutation: (i,j) is allowed
   iff i and the row matched to column j lie in one strong component of the switching digraph i -> row(j)).

   Mode V (values):      search_total_support V k H
     all u in {0..H}^k with u(0) = 0 (adding a constant changes nothing). Counts potentials with a tight permutation,
     with total support, and without; for those without total support, the histogram of the number of anti-diagonals
     carrying no allowed edge. Prints a line "COVER" for every potential with a tight permutation, without total
     support and with an allowed edge on every anti-diagonal (with its ray and facet flags), and a line "TSBIP" for
     every potential with total support whose tight graph has a bipartite component.
   Mode I (increments):  search_total_support I k H [LAST]
     all normalized potentials of height at most H: increments in {0..H}^(k-1) with minimum 0 and gcd 1, u(0) = 0
     (with LAST given, only those whose last increment equals LAST, to split a run). Counts the rays with a tight
     permutation and those without total support, and prints a line "FACET" for every facet ray without total
     support.

   Ranks: dim Z(E) = k - rank of the relations of Z, and the face dimension rank[B; X] - rank B over the allowed
   edges (B the row/column incidence matrix, totally unimodular, X the anti-diagonal rows). Ranks are computed modulo
   the primes 2^31 - 1 and 2^31 - 19 and the larger value is taken; a rank modulo a prime never exceeds the rational
   rank, so a reported ray (dim Z = 3) and a reported facet (face dimension 2k - 4) are certified. A ray or facet
   could be missed only if both modular ranks dropped below the rational rank.
   The last output line is a JSON summary.  Build: cc -O2 -o search_total_support search_total_support.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
#define MAXK 48
static int K;
static ll u[MAXK], v[2 * MAXK];
static unsigned char E[MAXK][MAXK], A[MAXK][MAXK];
static const ll PR1 = 2147483647LL, PR2 = 2147483629LL;

static ll powmod(ll a, ll e, ll p) { ll r = 1; a %= p; if (a < 0) a += p;
  while (e) { if (e & 1) r = (__int128)r * a % p; a = (__int128)a * a % p; e >>= 1; } return r; }

static int rank_mod(ll *M, int rows, int cols, ll p) {
  int r = 0;
  for (int c = 0; c < cols && r < rows; c++) {
    int piv = -1;
    for (int i = r; i < rows; i++) if (M[(ll)i * cols + c] % p) { piv = i; break; }
    if (piv < 0) continue;
    if (piv != r) for (int j = 0; j < cols; j++) { ll t = M[(ll)r * cols + j]; M[(ll)r * cols + j] = M[(ll)piv * cols + j]; M[(ll)piv * cols + j] = t; }
    ll inv = powmod(M[(ll)r * cols + c], p - 2, p);
    for (int j = 0; j < cols; j++) { ll x = M[(ll)r * cols + j] % p; if (x < 0) x += p; M[(ll)r * cols + j] = (__int128)x * inv % p; }
    for (int i = 0; i < rows; i++) if (i != r) {
      ll f = M[(ll)i * cols + c] % p; if (f < 0) f += p; if (!f) continue;
      for (int j = 0; j < cols; j++) { ll x = (M[(ll)i * cols + j] - (__int128)f * M[(ll)r * cols + j]) % p; if (x < 0) x += p; M[(ll)i * cols + j] = x; }
    }
    r++;
  }
  return r;
}

static int rank2(const ll *M, int rows, int cols) {
  ll *C = malloc(sizeof(ll) * (size_t)rows * cols + 8);
  memcpy(C, M, sizeof(ll) * (size_t)rows * cols); int a = rank_mod(C, rows, cols, PR1);
  memcpy(C, M, sizeof(ll) * (size_t)rows * cols); int b = rank_mod(C, rows, cols, PR2);
  free(C); return a > b ? a : b;
}

static void tight_graph(void) {
  for (int t = 0; t <= 2 * K - 2; t++) {
    ll m = (ll)9e18;
    for (int i = 0; i < K; i++) { int j = t - i; if (j < 0 || j >= K) continue; if (u[i] + u[j] < m) m = u[i] + u[j]; }
    v[t] = m;
  }
  for (int i = 0; i < K; i++) for (int j = 0; j < K; j++) E[i][j] = (u[i] + u[j] == v[i + j]);
}

/* perfect matching in E: mt[col] = row */
static int mt[MAXK], seen[MAXK];
static int augment(int i) {
  for (int j = 0; j < K; j++) if (E[i][j] && !seen[j]) { seen[j] = 1; if (mt[j] < 0 || augment(mt[j])) { mt[j] = i; return 1; } }
  return 0;
}
static int tight_permutation(void) {
  for (int j = 0; j < K; j++) mt[j] = -1;
  for (int i = 0; i < K; i++) { memset(seen, 0, sizeof seen); if (!augment(i)) return 0; }
  return 1;
}

/* strong components of the switching digraph i -> mt[j] for (i,j) in E (Tarjan) */
static int idx[MAXK], low[MAXK], onst[MAXK], stk[MAXK], sp, counter, comp[MAXK], ncomp;
static void tarjan(int x) {
  idx[x] = low[x] = counter++; stk[sp++] = x; onst[x] = 1;
  for (int j = 0; j < K; j++) if (E[x][j]) {
    int y = mt[j];
    if (idx[y] < 0) { tarjan(y); if (low[y] < low[x]) low[x] = low[y]; }
    else if (onst[y] && idx[y] < low[x]) low[x] = idx[y];
  }
  if (low[x] == idx[x]) { int y; do { y = stk[--sp]; onst[y] = 0; comp[y] = ncomp; } while (y != x); ncomp++; }
}
/* allowed graph; returns 1 iff total support */
static int allowed_graph(void) {
  for (int i = 0; i < K; i++) { idx[i] = -1; onst[i] = 0; }
  sp = counter = ncomp = 0;
  for (int i = 0; i < K; i++) if (idx[i] < 0) tarjan(i);
  int ts = 1;
  for (int i = 0; i < K; i++) for (int j = 0; j < K; j++) {
    A[i][j] = E[i][j] && comp[i] == comp[mt[j]];
    if (E[i][j] && !A[i][j]) ts = 0;
  }
  return ts;
}

static int empty_antidiagonals(void) {
  int empty = 0;
  for (int t = 0; t <= 2 * K - 2; t++) {
    int has = 0;
    for (int i = 0; i < K && !has; i++) { int j = t - i; if (j >= 0 && j < K && A[i][j]) has = 1; }
    empty += !has;
  }
  return empty;
}

/* dim Z(G): z(i) + z(j) constant over the pairs (i <= j) of G on each anti-diagonal */
static int dimZ(unsigned char G[MAXK][MAXK]) {
  int maxr = K * K; ll *M = calloc((size_t)maxr * K + 8, sizeof(ll)); int r = 0;
  for (int t = 0; t <= 2 * K - 2; t++) {
    int fi = -1, fj = -1;
    for (int i = 0; i < K; i++) {
      int j = t - i; if (j < i) break; if (j >= K || !G[i][j]) continue;
      if (fi < 0) { fi = i; fj = j; continue; }
      ll *row = M + (ll)r * K; row[i] += 1; row[j] += 1; row[fi] -= 1; row[fj] -= 1; r++;
    }
  }
  int rk = r ? rank2(M, r, K) : 0; free(M); return K - rk;
}

static int face_dim(void) {
  static int ei[MAXK * MAXK], ej[MAXK * MAXK]; int ne = 0;
  for (int i = 0; i < K; i++) for (int j = 0; j < K; j++) if (A[i][j]) { ei[ne] = i; ej[ne] = j; ne++; }
  int rows = 2 * K + 2 * K - 1; ll *M = calloc((size_t)rows * ne + 8, sizeof(ll));
  for (int e = 0; e < ne; e++) { M[(ll)ei[e] * ne + e] = 1; M[(ll)(K + ej[e]) * ne + e] = 1; M[(ll)(2 * K + ei[e] + ej[e]) * ne + e] = 1; }
  int rB = rank2(M, 2 * K, ne), rT = rank2(M, rows, ne); free(M); return rT - rB;
}

/* number of bipartite components (with at least one edge) of the graph of E */
static int bipartite_components_of_E(void) {
  int col[MAXK], b = 0;
  for (int i = 0; i < K; i++) col[i] = -1;
  for (int s = 0; s < K; s++) {
    if (col[s] >= 0) continue;
    int q[MAXK], h = 0, tl = 0, bip = 1; q[tl++] = s; col[s] = 0;
    while (h < tl) { int x = q[h++]; for (int y = 0; y < K; y++) if (E[x][y]) {
      if (col[y] < 0) { col[y] = 1 - col[x]; q[tl++] = y; } else if (col[y] == col[x]) bip = 0; } }
    b += bip;
  }
  return b;
}

static void print_u(const char *tag) { printf("%s", tag); for (int i = 0; i < K; i++) printf(" %lld", u[i]); }
static ll gcdl(ll a, ll b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { ll t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: search_total_support V|I k H [LAST]\n"); return 1; }
  char mode = argv[1][0]; K = atoi(argv[2]); int H = atoi(argv[3]); int last = argc >= 5 ? atoi(argv[4]) : -1;
  if (K < 3 || K > MAXK) { fprintf(stderr, "k out of range\n"); return 1; }
  if (mode == 'V') {
    long long total = 0, tp = 0, ts = 0, nots = 0, cover = 0, tsbip = 0, hist[2 * MAXK] = {0};
    int min_empty = 1 << 30;
    memset(u, 0, sizeof u);
    while (1) {
      total++; tight_graph();
      if (tight_permutation()) {
        tp++;
        if (allowed_graph()) {
          ts++;
          if (bipartite_components_of_E()) { tsbip++; print_u("TSBIP"); printf("\n"); }
        } else {
          nots++; int e = empty_antidiagonals(); hist[e]++; if (e < min_empty) min_empty = e;
          if (e == 0) { cover++; int ray = dimZ(E) == 3, fac = face_dim() == 2 * K - 4;
            print_u("COVER"); printf(" | ray %d facet %d\n", ray, fac); }
        }
      }
      int p = 1; while (p < K && u[p] == H) { u[p] = 0; p++; } if (p == K) break; u[p]++;
    }
    printf("{\"mode\": \"values\", \"k\": %d, \"H\": %d, \"potentials\": %lld, \"tight_permutation\": %lld, \"total_support\": %lld, "
           "\"no_total_support\": %lld, \"cover_without_total_support\": %lld, \"total_support_with_bipartite_component\": %lld, "
           "\"min_empty_antidiagonals\": %d, \"histogram_empty\": {", K, H, total, tp, ts, nots, cover, tsbip, nots ? min_empty : -1);
    int first = 1; for (int e = 0; e < 2 * MAXK; e++) if (hist[e]) { printf("%s\"%d\": %lld", first ? "" : ", ", e, hist[e]); first = 0; }
    printf("}}\n");
  } else if (mode == 'I') {
    int n = K - 1; int d[MAXK] = {0}; long long vectors = 0, rays = 0, rays_nots = 0, facets_nots = 0;
    while (1) {
      int mn = H + 1, mx = -1; ll g = 0;
      for (int i = 0; i < n; i++) { if (d[i] < mn) mn = d[i]; if (d[i] > mx) mx = d[i]; g = gcdl(g, d[i]); }
      if ((last < 0 || d[n - 1] == last) && mn == 0 && mx > 0 && g == 1) {
        vectors++; u[0] = 0; for (int i = 0; i < n; i++) u[i + 1] = u[i] + d[i];
        tight_graph();
        if (tight_permutation() && dimZ(E) == 3) {
          rays++;
          if (!allowed_graph()) { rays_nots++;
            if (empty_antidiagonals() == 0 && face_dim() == 2 * K - 4) { facets_nots++; print_u("FACET"); printf("\n"); fflush(stdout); } }
        }
      }
      int i = 0; while (i < n && d[i] == H) { d[i] = 0; i++; } if (i == n) break; d[i]++;
    }
    printf("{\"mode\": \"increments\", \"k\": %d, \"H\": %d, \"last\": %d, \"normalized_potentials\": %lld, \"rays_with_tight_permutation\": %lld, "
           "\"rays_without_total_support\": %lld, \"facet_rays_without_total_support\": %lld}\n", K, H, last, vectors, rays, rays_nots, facets_nots);
  } else { fprintf(stderr, "unknown mode\n"); return 1; }
  return 0;
}
