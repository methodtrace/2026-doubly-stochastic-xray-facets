/* Generator of the certificate trees for statement (S1) (checked by scripts/check_s1_certificates.py and
   scripts/check_s1_depth10.py; format in src/dsxray/s1_certificates.py). It extends attaining_structures.c by checked
   arithmetic, an optional time limit, export of every terminal record with its certificate, a fixed domain prefix,
   and connectivity pruning before the feasibility test.

   Usage for a certificate tree (positional arguments after those of attaining_structures.c):
     s1_certificate_tree S r -1 0 1 1 1 SECONDS PREFIX > tree.txt
   with SECONDS = 0 for no time limit and PREFIX = _ (the whole tree) or a comma-separated mask prefix such as
   1,1,3,2,7,5 (one domain of the r = 10 tree). Build: cc -O2 -o s1_certificate_tree s1_certificate_tree.c
   A stream is a proof only if it ends with its JSON summary and passes the Python checker. */
/* Exact enumeration of the attaining structures of integer (or real) sequences E_0, ..., E_r, and statement (S1).

   Usage: attaining_structures S r [split_level split_index split_count]   ((S1) mode)
          attaining_structures C r                                       (catalogue mode)
   Build: cc -O2 -o attaining_structures attaining_structures.c

   For a sequence E let M(s) = max_{p+q=s} E_p + E_q; the pair (p, q), p <= q, attains at level s = p + q if
   E_p + E_q = M(s). The attaining structure (AS) is the list, over the levels s = 0..2r-1, of the sets of attaining
   pairs (level 2r has the single pair (r, r)). The partners of q are the l with (q, l) attaining; q < r is solitary
   if it has at most one partner. The classes are the components of {0..r-1} under p ~ p' when p and p' both occur in
   attaining pairs of one level. (S1): if there are at least two classes, some q < r is solitary. (S1) and the
   classes depend on the AS only.

   The sequences with a given AS form a relatively open polyhedral cone: equalities between the attaining pairs of
   each level and strict inequalities "attaining > not attaining". The cone is invariant under E_p -> E_p + a p + b,
   which is removed by E_0 = E_1 = 0. A depth-first search over the levels chooses at each level a nonempty set of
   attaining pairs and keeps a branch only if the accumulated system is feasible. Feasibility is decided exactly:
   an integer basis of the nullspace of the equalities, then Gordan's alternative (y with C y > 0, or a nonnegative
   multiplier lambda != 0 with lambda C = 0) by a Phase-I simplex with Bland's rule in rational arithmetic on
   __int128 (the program aborts on overflow). Every "feasible" answer is backed by an integer witness E that is
   checked exactly against all constraints, every "infeasible" answer by a Gordan multiplier that is checked
   exactly. A witness found earlier is reused when it satisfies the new system.

   Modes
     attaining_structures S r [split_level split_index split_count]
        (S1) for all sequences with r boundary steps. Pruning: a branch is cut when some q < r can no longer reach
        two partners (no solitary index can then occur only if the leaf has every q < r with >= 2 partners), and
        when 0..r-1 already form a single class (classes only merge as levels are added). A leaf reached is a
        counterexample to (S1) and is printed as "S1 FAIL ...". With a split, only the branches at split_level whose
        running index is congruent to split_index modulo split_count are explored (to split a run into chunks).
     attaining_structures C r
        Catalogue: no pruning; every realizable AS is printed with its integer witness, as
        "AS 0:p,q;p,q;|1:...|...| W E_0 E_1 ... E_r".
   The last line is a JSON summary with the counts of nodes, feasibility tests, infeasibility certificates, leaves
   and, in catalogue mode, structures with several classes and those with every q < r having >= 2 partners.
   Build: cc -O2 -o attaining_structures attaining_structures.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

typedef __int128 i128;
#define MAXR 16
#define MAXROWS 512

static int R,CATALOGUE=0,early_connected=0,proof_records=0;
static void emit_zero(int row);
static int maskpath[2*MAXR],proof_depth=0;
static long long n_degree_pruned,n_rejected;
static int forced_mask[2*MAXR],forced_depth=0;
static void die(const char *m) { fprintf(stderr, "ABORT: %s\n", m); exit(2); }
static i128 negsafe(i128 a){i128 z;if(__builtin_sub_overflow((i128)0,a,&z))die("integer negation overflow");return z;}
static i128 iabs(i128 a) { return a < 0 ? negsafe(a) : a; }
static i128 gcd128(i128 a, i128 b) { a = iabs(a); b = iabs(b); while (b) { i128 t = a % b; a = b; b = t; } return a; }
static const i128 LIM = ((i128)1) << 100;

typedef struct { i128 n, d; } Q;
static void emit_gordan(int m,Q*lam);
static Q qn(i128 n, i128 d) {
  if (d == 0) die("division by zero");
  if (d < 0) { n = negsafe(n); d = negsafe(d); }
  i128 g = gcd128(n, d); if (g > 1) { n /= g; d /= g; }
  if (iabs(n) > LIM || d > LIM) die("overflow");
  Q q = {n, d}; return q;
}
static i128 mulsafe(i128 a,i128 b){i128 z;if(__builtin_mul_overflow(a,b,&z))die("integer multiplication overflow");return z;}
static i128 addsafe(i128 a,i128 b){i128 z;if(__builtin_add_overflow(a,b,&z))die("integer addition overflow");return z;}
static i128 subsafe(i128 a,i128 b){i128 z;if(__builtin_sub_overflow(a,b,&z))die("integer subtraction overflow");return z;}
static Q qadd(Q a,Q b){return qn(addsafe(mulsafe(a.n,b.d),mulsafe(b.n,a.d)),mulsafe(a.d,b.d));}
static Q qminus(Q a, Q b) { return qn(subsafe(mulsafe(a.n,b.d),mulsafe(b.n,a.d)),mulsafe(a.d,b.d)); }
static Q qmul(Q a, Q b) { return qn(mulsafe(a.n,b.n),mulsafe(a.d,b.d)); }
static Q qdiv(Q a, Q b) { return qn(mulsafe(a.n,b.d),mulsafe(a.d,b.n)); }
static int qsgn(Q a) { return a.n > 0 ? 1 : (a.n < 0 ? -1 : 0); }
static int qcmp(Q a, Q b) { i128 x=mulsafe(a.n,b.d),y=mulsafe(b.n,a.d); return x < y ? -1 : (x > y ? 1 : 0); }
static const Q QZERO = {0, 1}, QONE = {1, 1};

/* accumulated constraints on E (length R+1): EQ rows = 0, IN rows > 0 */
static int neq, nin;
static int EQ[MAXROWS][MAXR + 1], IN[MAXROWS][MAXR + 1];
static long long n_nodes, n_tests, n_certs, n_trivial, n_leaves, n_multi, n_alldeg2, n_conn_pruned;

/* integer basis (rows of B) of the nullspace of EQ plus the normalization E_0 = E_1 = 0 */
static int nullspace(i128 B[][MAXR + 1]) {
  int n = R + 1, m = neq + 2; static Q M[MAXROWS + 2][MAXR + 1];
  for (int i = 0; i < neq; i++) for (int j = 0; j < n; j++) M[i][j] = qn(EQ[i][j], 1);
  for (int j = 0; j < n; j++) { M[neq][j] = qn(j == 0, 1); M[neq + 1][j] = qn(j == 1, 1); }
  int piv[MAXR + 1], row = 0;
  for (int col = 0; col < n && row < m; col++) {
    int p = -1; for (int i = row; i < m; i++) if (M[i][col].n != 0) { p = i; break; }
    if (p < 0) continue;
    if (p != row) for (int j = 0; j < n; j++) { Q t = M[p][j]; M[p][j] = M[row][j]; M[row][j] = t; }
    Q pv = M[row][col]; for (int j = 0; j < n; j++) M[row][j] = qdiv(M[row][j], pv);
    for (int i = 0; i < m; i++) if (i != row && M[i][col].n != 0) {
      Q f = M[i][col]; for (int j = 0; j < n; j++) M[i][j] = qminus(M[i][j], qmul(f, M[row][j]));
    }
    piv[row] = col; row++;
  }
  int isp[MAXR + 1] = {0}; for (int i = 0; i < row; i++) isp[piv[i]] = 1;
  int k = 0;
  for (int fc = 0; fc < n; fc++) if (!isp[fc]) {
    Q v[MAXR + 1]; for (int j = 0; j < n; j++) v[j] = QZERO; v[fc] = QONE;
    for (int i = 0; i < row; i++) v[piv[i]] = qminus(QZERO, M[i][fc]);
    i128 L = 1; for (int j = 0; j < n; j++) { i128 g = gcd128(L, v[j].d); L=mulsafe(L/g,v[j].d); if (L > LIM) die("lcm overflow"); }
    for (int j = 0; j < n; j++) B[k][j] = mulsafe(v[j].n,L/v[j].d);
    k++;
  }
  return k;
}

/* Gordan: C is m x k (m inequality rows in nullspace coordinates). Returns 1 and an integer y with C y > 0 (checked),
   or 0 after checking a multiplier lambda >= 0, sum lambda = 1, lambda C = 0. */
static int gordan(int m, int k, i128 C[][MAXR + 1], i128 *yout) {
  int nrow = k + 1, ncol = m + nrow; static Q T[MAXR + 2][MAXROWS + MAXR + 3]; int basis[MAXR + 2];
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < m; j++) T[i][j] = qn(C[j][i], 1);
    for (int t = 0; t < nrow; t++) T[i][m + t] = qn(t == i, 1);
    T[i][ncol] = QZERO;
  }
  for (int j = 0; j < m; j++) T[k][j] = QONE;
  for (int t = 0; t < nrow; t++) T[k][m + t] = qn(t == k, 1);
  T[k][ncol] = QONE;
  for (int i = 0; i < nrow; i++) basis[i] = m + i;
  for (int iter = 0;; iter++) {     /* Phase I: minimize the sum of the artificial variables */
    if (iter > 100000) die("simplex iteration limit");
    int enter = -1;
    for (int j = 0; j < ncol && enter < 0; j++) {
      int inb = 0; for (int i = 0; i < nrow; i++) if (basis[i] == j) inb = 1;
      if (inb) continue;
      Q d = qn(j >= m, 1); for (int i = 0; i < nrow; i++) if (basis[i] >= m && T[i][j].n) d = qminus(d, T[i][j]);
      if (qsgn(d) < 0) enter = j;
    }
    if (enter < 0) break;
    int leave = -1; Q best = QZERO;
    for (int i = 0; i < nrow; i++) if (qsgn(T[i][enter]) > 0) {
      Q ratio = qdiv(T[i][ncol], T[i][enter]);
      if (leave < 0 || qcmp(ratio, best) < 0 || (qcmp(ratio, best) == 0 && basis[i] < basis[leave])) { best = ratio; leave = i; }
    }
    if (leave < 0) die("unbounded phase I");
    Q pv = T[leave][enter]; for (int j = 0; j <= ncol; j++) T[leave][j] = qdiv(T[leave][j], pv);
    for (int i = 0; i < nrow; i++) if (i != leave && T[i][enter].n) {
      Q f = T[i][enter]; for (int j = 0; j <= ncol; j++) T[i][j] = qminus(T[i][j], qmul(f, T[leave][j]));
    }
    basis[leave] = enter;
  }
  Q val = QZERO; for (int i = 0; i < nrow; i++) if (basis[i] >= m) val = qadd(val, T[i][ncol]);
  if (val.n == 0) {     /* lambda feasible: infeasibility of C y > 0; verify the multiplier */
    static Q lam[MAXROWS]; for (int j = 0; j < m; j++) lam[j] = QZERO;
    for (int i = 0; i < nrow; i++) if (basis[i] < m) lam[basis[i]] = T[i][ncol];
    Q s = QZERO;
    for (int j = 0; j < m; j++) { if (qsgn(lam[j]) < 0) die("negative multiplier"); s = qadd(s, lam[j]); }
    if (qcmp(s, QONE) != 0) die("multiplier sum");
    for (int i = 0; i < k; i++) {
      Q t = QZERO; for (int j = 0; j < m; j++) if (lam[j].n && C[j][i]) t = qadd(t, qmul(lam[j], qn(C[j][i], 1)));
      if (t.n) die("multiplier check");
    }
    n_certs++;emit_gordan(m,lam);return 0;
  }
  Q y[MAXR + 1];     /* dual values give y with C y > 0 */
  for (int t = 0; t < k; t++) {
    int j = m + t; Q d = QONE;
    for (int i = 0; i < nrow; i++) if (basis[i] >= m && T[i][j].n) d = qminus(d, T[i][j]);
    y[t] = qminus(QZERO, qminus(QONE, d));
  }
  i128 L = 1; for (int t = 0; t < k; t++) { i128 g = gcd128(L, y[t].d); L=mulsafe(L/g,y[t].d); if (L > LIM) die("lcm overflow"); }
  for (int t = 0; t < k; t++) yout[t] = mulsafe(y[t].n,L/y[t].d);
  for (int j = 0; j < m; j++) { i128 v = 0; for (int t = 0; t < k; t++) v=addsafe(v,mulsafe(C[j][t],yout[t])); if (v <= 0) die("witness not strict"); }
  return 1;
}

static i128 WIT[128][MAXR+1];static unsigned long long nwit=0; static i128 LASTW[MAXR + 1];
static int satisfies(const i128 *w) {
  for (int i = 0; i < neq; i++) { i128 v = 0; for (int j = 0; j <= R; j++) v=addsafe(v,mulsafe(EQ[i][j],w[j])); if (v) return 0; }
  for (int i = 0; i < nin; i++) { i128 v = 0; for (int j = 0; j <= R; j++) v=addsafe(v,mulsafe(IN[i][j],w[j])); if (v <= 0) return 0; }
  return w[0] == 0 && w[1] == 0;
}
static int feasible(void) {
  for(unsigned t=0;t<64&&t<nwit;t++){unsigned h=(nwit-1-t)%128;if(satisfies(WIT[h])){memcpy(LASTW,WIT[h],sizeof LASTW);return 1;}}
  n_tests++;
  static i128 B[MAXR + 1][MAXR + 1]; int k = nullspace(B);
  if (nin == 0) { memset(LASTW, 0, sizeof LASTW); if (!satisfies(LASTW)) die("zero witness"); return 1; }
  if(k==0){emit_zero(0);return 0;}
  static i128 C[MAXROWS][MAXR + 1];
  for (int j = 0; j < nin; j++) {
    int zero = 1;
    for (int t = 0; t < k; t++) { i128 v = 0; for (int c = 0; c <= R; c++) v=addsafe(v,mulsafe(IN[j][c],B[t][c])); C[j][t] = v; if (v) zero = 0; }
    if(zero){n_trivial++;emit_zero(j);return 0;}   /* an inequality row vanishing on the nullspace reads 0 > 0 */
  }
  i128 y[MAXR + 1]; if (!gordan(nin, k, C, y)) return 0;
  i128 x[MAXR + 1];
  for (int c = 0; c <= R; c++) { x[c] = 0; for (int t = 0; t < k; t++) x[c]=addsafe(x[c],mulsafe(y[t],B[t][c])); if (iabs(x[c]) > LIM) die("witness overflow"); }
  if (!satisfies(x)) die("witness check failed");
  memcpy(WIT[nwit%128],x,sizeof x);nwit++;
  memcpy(LASTW, x, sizeof LASTW);
  return 1;
}

static void print128(i128 x){if(x<0){putchar('-');x=negsafe(x);}char b[64];int n=0;do{b[n++]='0'+x%10;x/=10;}while(x);while(n)putchar(b[--n]);}
static void print_prefix(char kind){printf("%c ",kind);for(int s=0;s<proof_depth;s++)printf("%s%d",s?",":"",maskpath[s]);if(!proof_depth)putchar('_');}
static void emit_zero(int row){n_rejected++;if(proof_records){print_prefix('Z');printf(" %d\n",row);}}
static void emit_gordan(int m,Q*lam){n_rejected++;if(proof_records){print_prefix('G');for(int j=0;j<m;j++)if(lam[j].n){printf(" %d:",j);print128(lam[j].n);putchar('/');print128(lam[j].d);}putchar('\n');}}
static void emit_degree(int q){n_degree_pruned++;if(proof_records){print_prefix('D');printf(" %d\n",q);}}
static void emit_connected(void){if(proof_records){print_prefix('C');putchar('\n');}}
static int deg[MAXR + 1], remaining[MAXR + 1], mindeg = 2, prune_conn = 1;
static int chosen[2 * MAXR][MAXR + 1][2], nch[2 * MAXR];
static int split_level = -1, split_index = 0, split_count = 1; static long long split_counter = 0;
static int pairs_at(int s, int P[][2]) { int n = 0; for (int a = (s - R > 0 ? s - R : 0); a <= s / 2; a++) { P[n][0] = a; P[n][1] = s - a; n++; } return n; }
static int par[MAXR + 1];
static int find(int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; }
static int classes_upto(int smax) {   /* number of classes of 0..R-1 from the levels 0..smax */
  for (int p = 0; p < R; p++) par[p] = p;
  for (int s = 0; s <= smax; s++) {
    int first = -1;
    for (int t = 0; t < nch[s]; t++) for (int e = 0; e < 2; e++) {
      int p = chosen[s][t][e];
      if (p < R) { if (first < 0) first = p; else { int a = find(p), b = find(first); if (a != b) par[a] = b; } }
    }
  }
  int c = 0; for (int p = 0; p < R; p++) c += find(p) == p; return c;
}
static void print_structure(void) {
  printf("AS ");
  for (int s = 0; s < 2 * R; s++) { printf("%d:", s); for (int t = 0; t < nch[s]; t++) printf("%d,%d;", chosen[s][t][0], chosen[s][t][1]); printf("|"); }
  printf(" W"); for (int j = 0; j <= R; j++) {putchar(' ');print128(LASTW[j]);} printf("\n");
}
static void leaf(void) {
  n_leaves++;
  int multi = classes_upto(2 * R - 1) > 1, all2 = 1;
  for (int q = 0; q < R; q++) if (deg[q] < 2) all2 = 0;
  if (CATALOGUE) { print_structure(); n_multi += multi; n_alldeg2 += all2; if (multi && all2) printf("S1 FAIL\n"); return; }
  if (multi && all2) {
    n_multi++;
    if(proof_records){proof_depth=2*R;print_prefix('X');for(int j=0;j<=R;j++){putchar(' ');print128(LASTW[j]);}putchar('\n');return;}
    printf("S1 FAIL: structure");
    for (int s = 0; s < 2 * R; s++) { printf(" [%d:", s); for (int t = 0; t < nch[s]; t++) printf("(%d,%d)", chosen[s][t][0], chosen[s][t][1]); printf("]"); }
    printf("\n"); fflush(stdout);
  }
}
static void addrow(int (*M)[MAXR + 1], int *n, int a, int b, int c, int d) {
  if (*n >= MAXROWS) die("too many rows");
  memset(M[*n], 0, sizeof(int) * (MAXR + 1)); M[*n][a]++; M[*n][b]++; M[*n][c]--; M[*n][d]--; (*n)++;
}
static void dfs(int s) {
  n_nodes++;
  if(n_nodes%100000==0)fprintf(stderr,"progress r=%d nodes=%lld tests=%lld certs=%lld cpu_seconds=%.1f\n",R,n_nodes,n_tests,n_certs,(double)clock()/CLOCKS_PER_SEC);
  if (s == 2 * R) { leaf(); return; }
  int P[MAXR + 1][2]; int n = pairs_at(s, P);
  int parent_count=classes_upto(s-1),parent_save[MAXR+1];memcpy(parent_save,par,sizeof par);
  for (int mask = 1; mask < (1 << n); mask++) {
    if(s<forced_depth && mask!=forced_mask[s])continue;
    maskpath[s]=mask;proof_depth=s+1;
    if (s == split_level) { long long id = split_counter++; if (id % split_count != split_index) continue; }
    for (int i = 0; i < n; i++) { remaining[P[i][0]]--; if (P[i][1] != P[i][0]) remaining[P[i][1]]--; }
    for (int i = 0; i < n; i++) if (mask >> i & 1) { deg[P[i][0]]++; if (P[i][1] != P[i][0]) deg[P[i][1]]++; }
    int ok=1;for(int q=0;q<R;q++)if(deg[q]+remaining[q]<mindeg){ok=0;emit_degree(q);break;}
    if (ok) {
      int f = -1; for (int i = 0; i < n; i++) if (mask >> i & 1) { f = i; break; }
      nch[s]=0;
      for(int i=0;i<n;i++)if(mask>>i&1){chosen[s][nch[s]][0]=P[i][0];chosen[s][nch[s]][1]=P[i][1];nch[s]++;}
      memcpy(par,parent_save,sizeof par);int current_classes=parent_count,first_index=-1;
      for(int t=0;t<nch[s];t++)for(int z=0;z<2;z++){
        int p=chosen[s][t][z];if(p>=R)continue;
        if(first_index<0)first_index=p;else{int a=find(p),b=find(first_index);if(a!=b){par[a]=b;current_classes--;}}
      }
      int earlyprune=early_connected&&prune_conn&&current_classes==1;
      if(earlyprune){n_conn_pruned++;emit_connected();}
      int e0 = neq, i0 = nin;
      for (int i = 0; i < n && !earlyprune; i++) if (i != f) {
        if (mask >> i & 1) addrow(EQ, &neq, P[f][0], P[f][1], P[i][0], P[i][1]);
        else addrow(IN, &nin, P[f][0], P[f][1], P[i][0], P[i][1]);
      }
      if (!earlyprune && feasible()) {
        nch[s] = 0;
        for (int i = 0; i < n; i++) if (mask >> i & 1) { chosen[s][nch[s]][0] = P[i][0]; chosen[s][nch[s]][1] = P[i][1]; nch[s]++; }
        if(prune_conn&&current_classes==1){n_conn_pruned++;emit_connected();}
        else dfs(s + 1);
      }
      neq = e0; nin = i0;
    }
    for (int i = 0; i < n; i++) if (mask >> i & 1) { deg[P[i][0]]--; if (P[i][1] != P[i][0]) deg[P[i][1]]--; }
    for (int i = 0; i < n; i++) { remaining[P[i][0]]++; if (P[i][1] != P[i][0]) remaining[P[i][1]]++; }
  }
}

static void timed_out(int sig){(void)sig;ssize_t z=write(2,"TIME_LIMIT\n",11);(void)z;_exit(124);}
int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: attaining_structures S|C r [split_level split_index split_count]\n"); return 1; }
  CATALOGUE = argv[1][0] == 'C'; R = atoi(argv[2]);
  if (R < 1 || R >= MAXR) die("r out of range");
  if (CATALOGUE) { mindeg = 0; prune_conn = 0; }
  if(argc>=6){split_level=atoi(argv[3]);split_index=atoi(argv[4]);split_count=atoi(argv[5]);}
  if(split_count<1||split_index<0||split_index>=split_count||split_level<-1||split_level>=2*R)die("invalid legacy split arguments");
  if(argc>=7)proof_records=atoi(argv[6]);
  if(argc>=8)early_connected=atoi(argv[7]);
  if(argc>=9){int seconds=atoi(argv[8]);if(seconds>0){signal(SIGALRM,timed_out);alarm(seconds);}}
  if(argc>=10 && strcmp(argv[9],"_")){
    char buf[512];if(strlen(argv[9])>=sizeof buf)die("domainprefix too long");strcpy(buf,argv[9]);
    for(char*t=strtok(buf,",");t;t=strtok(NULL,",")){if(forced_depth>=2*R)die("domainprefix too deep");forced_mask[forced_depth++]=atoi(t);}
  }
  for(int s=0;s<forced_depth;s++){int P[MAXR+1][2],n=pairs_at(s,P);if(forced_mask[s]<1||forced_mask[s]>=(1<<n))die("invalid domainprefix mask");}
  if(proof_records&&(CATALOGUE||split_count!=1))die("proof export requires S1mode and fixed-prefix splitting");
  if(proof_records){printf("T {\"format\":\"s1-tree-v1\",\"r\":%d,\"early_connected\":%d,\"domain_prefix\":[",R,early_connected);for(int j=0;j<forced_depth;j++)printf("%s%d",j?",":"",forced_mask[j]);printf("]}\n");}
  for (int s = 0; s < 2 * R; s++) {
    int P[MAXR + 1][2]; int n = pairs_at(s, P);
    for (int i = 0; i < n; i++) { remaining[P[i][0]]++; if (P[i][1] != P[i][0]) remaining[P[i][1]]++; }
  }
  dfs(0);
  if(proof_records){printf("{\"mode\":\"s1-proof\",\"r\":%d,\"early_connected\":%d,\"nodes\":%lld,\"feasibility_tests\":%lld,\"infeasibility_certificates\":%lld,\"trivial_infeasibilities\":%lld,\"rejected_branches\":%lld,\"degree_pruned\":%lld,\"connectivity_pruned\":%lld,\"leaves\":%lld,\"s1_failures\":%lld}\n",R,early_connected,n_nodes,n_tests,n_certs,n_trivial,n_rejected,n_degree_pruned,n_conn_pruned,n_leaves,n_multi);return n_multi?3:0;}
  printf("{\"mode\": \"%s\", \"r\": %d, \"split\": [%d, %d, %d], \"nodes\": %lld, \"feasibility_tests\": %lld, "
         "\"infeasibility_certificates\": %lld, \"trivial_infeasibilities\": %lld, \"rejected_branches\": %lld, \"degree_pruned\": %lld, \"connectivity_pruned\": %lld, \"leaves\": %lld, \"several_classes\": %lld, "
         "\"all_degrees_at_least_two\": %lld, \"s1_failures\": %lld}\n",
         CATALOGUE ? "catalogue" : "s1", R, split_level, split_index, split_count, n_nodes, n_tests, n_certs,
         n_trivial, n_rejected, n_degree_pruned, n_conn_pruned, n_leaves, n_multi, n_alldeg2, CATALOGUE ? 0LL : n_multi);
  return 0;
}
