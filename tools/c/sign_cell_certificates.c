/* Generator of the sign-cell certificate streams for statement E-flat (checked by scripts/check_cells_depth7.py;
   format in src/dsxray/sign_certificates.py).

   Usage: sign_cell_certificates r 1 0 1 > stream.txt      (records on, empty unit set only: the two-cell reduction)
   Build: cc -O2 -o sign_cell_certificates sign_cell_certificates.c   (needs __int128 and the overflow builtins)
   A stream is a proof only if it ends with its JSON summary and passes the Python checker. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
#define MAXR 16
#define MAXROWS 512

static int R;
static void emit_zero(int row);

static void die(const char *m) { fprintf(stderr, "ABORT: %s\n", m); exit(2); }
static i128 negsafe(i128 a){i128 z;if(__builtin_sub_overflow((i128)0,a,&z))die("integer negation overflow");return z;}
static i128 iabs(i128 a) { return a < 0 ? negsafe(a) : a; }
static i128 gcd128(i128 a, i128 b) { a = iabs(a); b = iabs(b); while (b) { i128 t = a % b; a = b; b = t; } return a; }
static const i128 LIM = ((i128)1) << 100;

typedef struct { i128 n, d; } Q;
static void emit_gordan(int m,Q* lam);
static Q qn(i128 n, i128 d) {
  if (d == 0) die("division by zero");
  if (d < 0) { n = negsafe(n); d = negsafe(d); }
  i128 g = gcd128(n, d); if (g > 1) { n /= g; d /= g; }
  if (iabs(n) > LIM || d > LIM) die("overflow");
  Q q = {n, d}; return q;
}
static i128 mulsafe(i128 a,i128 b) {i128 z;if(__builtin_mul_overflow(a,b,&z))die("integer multiplication overflow");return z;}
static i128 addsafe(i128 a,i128 b) {i128 z;if(__builtin_add_overflow(a,b,&z))die("integer addition overflow");return z;}
static i128 subsafe(i128 a,i128 b) {i128 z;if(__builtin_sub_overflow(a,b,&z))die("integer subtraction overflow");return z;}
static Q qadd(Q a, Q b) { return qn(addsafe(mulsafe(a.n,b.d),mulsafe(b.n,a.d)),mulsafe(a.d,b.d)); }
static Q qminus(Q a, Q b) { return qn(subsafe(mulsafe(a.n,b.d),mulsafe(b.n,a.d)),mulsafe(a.d,b.d)); }
static Q qmul(Q a, Q b) { return qn(mulsafe(a.n,b.n),mulsafe(a.d,b.d)); }
static Q qdiv(Q a, Q b) { return qn(mulsafe(a.n,b.d),mulsafe(a.d,b.n)); }
static int qsgn(Q a) { return a.n > 0 ? 1 : (a.n < 0 ? -1 : 0); }
static int qcmp(Q a, Q b) { i128 x = mulsafe(a.n,b.d), y = mulsafe(b.n,a.d); return x < y ? -1 : (x > y ? 1 : 0); }
static const Q QZERO = {0, 1}, QONE = {1, 1};

/* accumulated constraints on E (length R+1): EQ rows = 0, IN rows > 0 */
static int neq, nin;
static int EQ[MAXROWS][MAXR + 1], IN[MAXROWS][MAXR + 1];
static long long n_nodes, n_tests, n_certs, n_trivial, n_leaves;

/* integer basis (rows of B) of the nullspace of EQ plus the normalization E_0 = E_1 = 0 */
static int nullspace(i128 B[][MAXR + 1]) {
  int n = R, m = neq; static Q M[MAXROWS + 2][MAXR + 1];
  for (int i = 0; i < neq; i++) for (int j = 0; j < n; j++) M[i][j] = qn(EQ[i][j], 1);

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
    i128 L = 1; for (int j = 0; j < n; j++) { i128 g = gcd128(L, v[j].d); L = mulsafe(L / g, v[j].d); if (L > LIM) die("lcm overflow"); }
    for (int j = 0; j < n; j++) B[k][j] = mulsafe(v[j].n, (L / v[j].d));
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
    n_certs++; emit_gordan(m,lam); return 0;
  }
  Q y[MAXR + 1];     /* dual values give y with C y > 0 */
  for (int t = 0; t < k; t++) {
    int j = m + t; Q d = QONE;
    for (int i = 0; i < nrow; i++) if (basis[i] >= m && T[i][j].n) d = qminus(d, T[i][j]);
    y[t] = qminus(QZERO, qminus(QONE, d));
  }
  i128 L = 1; for (int t = 0; t < k; t++) { i128 g = gcd128(L, y[t].d); L = mulsafe(L / g, y[t].d); if (L > LIM) die("lcm overflow"); }
  for (int t = 0; t < k; t++) yout[t] = mulsafe(y[t].n, (L / y[t].d));
  for (int j = 0; j < m; j++) { i128 v = 0; for (int t = 0; t < k; t++) v=addsafe(v,mulsafe(C[j][t],yout[t])); if (v <= 0) die("witness not strict"); }
  return 1;
}

static i128 WIT[128][MAXR + 1]; static unsigned long long nwit = 0; static i128 LASTW[MAXR + 1];
static int satisfies(const i128 *w) {
  for (int i = 0; i < neq; i++) { i128 v = 0; for (int j = 0; j < R; j++) v=addsafe(v,mulsafe(EQ[i][j],w[j])); if (v) return 0; }
  for (int i = 0; i < nin; i++) { i128 v = 0; for (int j = 0; j < R; j++) v=addsafe(v,mulsafe(IN[i][j],w[j])); if (v <= 0) return 0; }
  return 1;
}
static int feasible(void) {
  for (unsigned t=0;t<64 && t<nwit;t++) {unsigned h=(nwit-1-t)%128;if(satisfies(WIT[h])){memcpy(LASTW,WIT[h],sizeof LASTW);return 1;}}
  n_tests++;
  static i128 B[MAXR + 1][MAXR + 1]; int k = nullspace(B);
  if (nin == 0) { memset(LASTW, 0, sizeof LASTW); if (!satisfies(LASTW)) die("zero witness"); return 1; }
  if (k == 0) {emit_zero(0);return 0;}
  static i128 C[MAXROWS][MAXR + 1];
  for (int j = 0; j < nin; j++) {
    int zero = 1;
    for (int t = 0; t < k; t++) { i128 v = 0; for (int c = 0; c < R; c++) v=addsafe(v,mulsafe(IN[j][c],B[t][c])); C[j][t] = v; if (v) zero = 0; }
    if (zero) { n_trivial++; emit_zero(j);return 0; }   /* an inequality row vanishing on the nullspace reads 0 > 0 */
  }
  i128 y[MAXR + 1]; if (!gordan(nin, k, C, y)) return 0;
  i128 x[MAXR + 1];
  for (int c = 0; c < R; c++) { x[c] = 0; for (int t = 0; t < k; t++) x[c]=addsafe(x[c],mulsafe(y[t],B[t][c])); if (iabs(x[c]) > LIM) die("witness overflow"); }
  if (!satisfies(x)) die("witness check failed");
  memcpy(WIT[nwit%128],x,sizeof x);nwit++;
  memcpy(LASTW, x, sizeof LASTW);
  return 1;
}


/* Sign-cell enumeration for E-flat. All comparison forms have coefficient sum
   zero. Thus x_i=1 on I, x_i>=2 otherwise is equivalent, up to scaling x-1,
   to homogeneous y_i=0 on I, y_i>0 otherwise. */
static int forms[256][MAXR+1], nf, sgpath[256], maskI, write_records;
static long long counts[4], reject_count;
static void print128(i128 x){
  if(x<0){putchar('-');x=negsafe(x);} char buf[64];int n=0;
  do {buf[n++]='0'+x%10;x/=10;}while(x);
  while(n)putchar(buf[--n]);
}
static void prefix_record(char kind,int depth){
  printf("%c %d ",kind,maskI);
  for(int i=0;i<depth;i++)putchar(sgpath[i]>0?'+':sgpath[i]<0?'-':'0');
  if(depth==0)putchar('_');
}
static int cert_depth;
static void emit_zero(int row){
  reject_count++;
  if(write_records){prefix_record('F',cert_depth);printf(" Z %d\n",row);}
}
static void emit_gordan(int m,Q*lam){
  reject_count++;
  if(write_records){prefix_record('F',cert_depth);printf(" G");
    for(int j=0;j<m;j++)if(lam[j].n){printf(" %d:",j);print128(lam[j].n);putchar('/');print128(lam[j].d);}
    putchar('\n');}
}
static void gen_forms(void){
  for(int s=0;s<=2*R;s++) {
    int P[MAXR+1][2],n=0;
    for(int a=(s-R>0?s-R:0);a<=s/2;a++){P[n][0]=a;P[n][1]=s-a;n++;}
    for(int p=0;p<n;p++)for(int q=p+1;q<n;q++){
      int f[MAXR+1]={0},orient=0;
      for(int j=0;j<R;j++)f[j]=(j<P[p][0])+(j<P[p][1])-(j<P[q][0])-(j<P[q][1]);
      for(int j=0;j<R;j++)if(f[j]){orient=f[j]>0?1:-1;break;}
      if(!orient)die("zero comparisonform");
      int total=0;for(int j=0;j<R;j++){f[j]*=orient;total+=f[j];}
      if(total)die("comparisonform coefficients do notsumzero");
      int old=0;for(int h=0;h<nf;h++)if(!memcmp(f,forms[h],R*sizeof(int)))old=1;
      if(!old){if(nf>=256)die("too manyforms");memcpy(forms[nf++],f,sizeof f);}
    }
  }
}
static unsigned adj[MAXR],starK,starA,starB;
static int matched[MAXR+1];
static int aug(int q,unsigned*seen){
  for(int c=0;c<=R;c++)if((adj[q]>>c&1)&&!(*seen>>c&1)){
    *seen|=1u<<c;
    if(matched[c]<0||aug(matched[c],seen)){matched[c]=q;return 1;}
  }return 0;
}
static int verdict(void){
  i128 E[MAXR+1],M[2*MAXR+1];int nmax[2*MAXR+1];
  E[0]=0;for(int i=0;i<R;i++)E[i+1]=addsafe(E[i],LASTW[i]);
  for(int s=0;s<=2*R;s++){
    M[s]=-1;nmax[s]=0;
    for(int p=(s-R>0?s-R:0);p<=s/2;p++){
      i128 z=addsafe(E[p],E[s-p]);
      if(z>M[s]){M[s]=z;nmax[s]=1;}else if(z==M[s])nmax[s]++;
    }
  }
  for(int q=0;q<R;q++){
    adj[q]=0;for(int c=0;c<=R;c++)if(addsafe(E[q],E[c])==M[q+c])adj[q]|=1u<<c;
  }
  for(int c=0;c<=R;c++)matched[c]=-1;
  for(int q=0;q<R;q++){unsigned seen=0;if(!aug(q,&seen))return 0;}
  unsigned reach=0;int changed=1;
  while(changed){changed=0;for(int q=0;q<R;q++)if(!(reach>>q&1)){
    for(int c=0;c<=R;c++)if((adj[q]>>c&1)&&(matched[c]<0||(reach>>matched[c]&1))){reach|=1u<<q;changed=1;break;}
  }}
  starK=((1u<<R)-1)^reach;if(!starK)return 1;
  starA=0;for(int q=0;q<R;q++)if(starK>>q&1)starA|=adj[q];
  if(__builtin_popcount(starA)!=__builtin_popcount(starK))die("balancedsetcount");
  starB=0;for(int d=1;d<=R;d++)if((maskI>>(d-1)&1)&&(starK>>(d-1)&1))starB|=1u<<d;
  for(int c=0;c<=R;c++)if(starA>>c&1)for(int d=0;d<=R;d++){
    if(!(starB>>d&1)&&nmax[c+d]==1&&addsafe(E[c],E[d])==M[c+d])return 2;
  }
  return 3;
}
static void cell_leaf(void){
  n_leaves++;int v=verdict();counts[v]++;
  if(write_records||v==3){prefix_record('L',nf);printf(" %d",v);for(int j=0;j<R;j++){putchar(' ');print128(LASTW[j]);}putchar('\n');}
}
static void cells_dfs(int depth){
  n_nodes++;
  if(depth==nf){cell_leaf();return;}
  int saved_eq=neq,saved_in=nin;
  const int signs[3]={1,-1,0};
  for(int b=0;b<3;b++){
    int s=signs[b];sgpath[depth]=s;cert_depth=depth+1;
    int(*mat)[MAXR+1]=s?IN:EQ;int*nr=s?&nin:&neq;
    if(*nr>=MAXROWS)die("too many constraints");
    for(int j=0;j<R;j++)mat[*nr][j]=forms[depth][j]*(s?s:1);
    (*nr)++;
    if(feasible())cells_dfs(depth+1);
    neq=saved_eq;nin=saved_in;
  }
}
int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: cells_exact r [records] [unitmask_start unitmask_stop]\n");return 1;}
  R=atoi(argv[1]);write_records=argc>=3?atoi(argv[2]):0;
  if(R<1||R>=MAXR)die("r outofrange");
  gen_forms();int start=argc>=4?atoi(argv[3]):0,stop=argc>=5?atoi(argv[4]):1<<R;
  for(maskI=start;maskI<stop;maskI++){
    neq=nin=0;
    for(int j=0;j<R;j++){
      int(*mat)[MAXR+1]=(maskI>>j&1)?EQ:IN;int*nr=(maskI>>j&1)?&neq:&nin;
      memset(mat[*nr],0,(MAXR+1)*sizeof(int));mat[*nr][j]=1;(*nr)++;
    }
    cert_depth=0;if(!feasible())die("coordinatepositivebaseinfeasible");
    cells_dfs(0);
    fprintf(stderr,"progress r=%d mask=%d cells=%lld tests=%lld rejects=%lld\n",R,maskI,n_leaves,n_tests,reject_count);
  }
  printf("{\"r\":%d,\"forms\":%d,\"unitmask_start\":%d,\"unitmask_stop\":%d,\"nodes\":%lld,\"feasibility_tests\":%lld,\"rejected_branches\":%lld,\"cells\":%lld,\"hall_fails\":%lld,\"kstar_empty\":%lld,\"eflat_holds\":%lld,\"eflat_fails\":%lld}\n",R,nf,start,stop,n_nodes,n_tests,reject_count,n_leaves,counts[0],counts[1],counts[2],counts[3]);
  return counts[3]?3:0;
}
