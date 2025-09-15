// red.c — Recombination-aware Edit Distance (RED)
// Author: Katsuhiko Hayashi
// Build:   gcc -O3 -std=c11 -o red red.c
//
// Usage:
//   ./red --parents-aln parents.aligned.fasta --child child.fasta \
//           --out out/prefix \
//           [--mismatch 1] [--gap 2] [--switch 1] \
//           [--win 21] [--kappa 1] [--stick 0] \
//           [--mu-inf <int>] [--mu-eq <int>] \
//           [--gap-swell 1] [--gap-near 10] [--reward 0]
//
// Notes:
// - Parents: exactly two aligned sequences (equal length), may contain '-'.
// - Child:   one or more ungapped sequences.
// - Outputs per child:
//     <outprefix>.<childName>.origins.tsv   (j \t parent[A|B|INS])
//     <outprefix>.<childName>.summary.tsv   (total,mismatches,insertions,deletions,switches)
//
// Model (brief):
// - Informative site weighting: mismatch cost μ_t is higher when A[t]!=B[t] (non-gap).
// - Gap-aware switch: λ'(t)=λ_open + δ_gap if near gaps within ±g in parents.
// - Stickiness reward: staying in same parent on a diagonal step adds -σ.
// - Bootscan unary bias: β_A(j)=-κ s(j), β_B(j)=+κ s(j), where s(j) = matches_A - matches_B
//   over a window of size W on ungapped parents, counting only informative positions.
// - β_m(j) is added once on transitions that consume the child (diag/switch/ins), not on deletions.

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <limits.h>

typedef struct { char *name; char *seq; size_t len; } FastaRec;
typedef struct { FastaRec *recs; size_t n, cap; } FastaVec;

static void die(const char *msg){ fprintf(stderr,"[ERROR] %s\n", msg); exit(1); }
static void *xmalloc(size_t n){ void *p=malloc(n); if(!p) die("malloc failed"); return p; }
static void *xcalloc(size_t n, size_t sz){ void *p=calloc(n,sz); if(!p) die("calloc failed"); return p; }
static void *xrealloc(void *q, size_t n){ void *p=realloc(q,n); if(!p) die("realloc failed"); return p; }
static char *xstrdup(const char *s){ char *p=strdup(s); if(!p) die("strdup failed"); return p; }

static void to_upper_inplace(char *s){ for(;*s;++s)*s=(char)toupper((unsigned char)*s); }

static void fasta_vec_init(FastaVec *v){ v->recs=NULL; v->n=0; v->cap=0; }
static void fasta_vec_push(FastaVec *v, FastaRec r){
  if(v->n==v->cap){ v->cap = v->cap? v->cap*2 : 8; v->recs = xrealloc(v->recs, v->cap*sizeof(FastaRec)); }
  v->recs[v->n++] = r;
}

// Read multi-FASTA (concatenate wrapped lines; strip whitespace)
static void read_fasta(const char *path, FastaVec *out){
  fasta_vec_init(out);
  FILE *fp=fopen(path,"r"); if(!fp){ perror(path); die("open fasta failed"); }
  char *line=NULL; size_t cap=0; ssize_t n; char *name=NULL;
  char *buf=NULL; size_t blen=0, bcap=0;
  while((n=getline(&line,&cap,fp))!=-1){
    if(n>0 && (line[n-1]=='\n'||line[n-1]=='\r')) line[--n]=0;
    if(!n) continue;
    if(line[0]=='>'){
      if(name){
        FastaRec rec={ xstrdup(name), xstrdup(buf?buf:""), blen };
        fasta_vec_push(out, rec);
        blen=0; if(buf) buf[0]=0;
      }
      name=line+1;
    }else{
      if(!buf){ bcap=(size_t)n+64; buf=xmalloc(bcap); buf[0]=0; }
      if(blen+(size_t)n+1>bcap){ bcap=(blen+(size_t)n+1)*2; buf=xrealloc(buf,bcap); }
      for(size_t i=0;i<(size_t)n;i++){ if(!isspace((unsigned char)line[i])) buf[blen++]=line[i]; }
      buf[blen]=0;
    }
  }
  if(name){ FastaRec rec={ xstrdup(name), xstrdup(buf?buf:""), blen }; fasta_vec_push(out, rec); }
  free(line); if(buf) free(buf); fclose(fp);
  for(size_t i=0;i<out->n;i++){ to_upper_inplace(out->recs[i].seq); out->recs[i].len=strlen(out->recs[i].seq); }
}

static char *ungap(const char *s, size_t *out_len){
  size_t n=strlen(s);
  char *u=xmalloc(n+1); size_t k=0;
  for(size_t i=0;i<n;i++){ if(s[i]!='-') u[k++]=s[i]; }
  u[k]=0; if(out_len) *out_len=k; return u;
}

// ---------- CLI ----------
static const char *get_arg(int argc, char **argv, const char *key){
  for(int i=1;i<argc;i++){ if(strcmp(argv[i],key)==0){ if(i+1<argc && argv[i+1][0]!='-') return argv[i+1]; else return "true"; } }
  return NULL;
}

static void usage(){
  fprintf(stderr,
    "Usage: red --parents-aln PARENTS.fasta --child CHILD.fasta --out OUTPREFIX\n"
    "             [--mismatch 1] [--gap 2] [--switch 1]\n"
    "             [--win 21] [--kappa 1] [--stick 0]\n"
    "             [--mu-inf <int>] [--mu-eq <int>]\n"
    "             [--gap-swell 1] [--gap-near 10] [--reward 0]\n");
}

// ---------- Output helpers ----------
static void write_origins_tsv(const char *path, const char *orig, size_t J){
  FILE *fp=fopen(path,"w"); if(!fp){ perror(path); die("open origins out failed"); }
  fprintf(fp,"j\tparent\n");
  for(size_t j=1;j<=J;j++){
    char c=orig[j-1];
    if(c=='I') fprintf(fp,"%zu\tINS\n", j);
    else       fprintf(fp,"%zu\t%c\n",  j, c); // 'A' or 'B'
  }
  fclose(fp);
}
static void write_summary_tsv(const char *path, int total, int mism, int ins, int del, int sw){
  FILE *fp=fopen(path,"w"); if(!fp){ perror(path); die("open summary out failed"); }
  fprintf(fp,"total\tmismatches\tinsertions\tdeletions\tswitches\n");
  fprintf(fp,"%d\t%d\t%d\t%d\t%d\n", total, mism, ins, del, sw);
  fclose(fp);
}

// ---------- Backpointers ----------
enum { OP_NONE=0, OP_DIAG=1, OP_INS=2, OP_DEL=3, OP_SWITCH=4 };
typedef struct {
  uint8_t *opA, *pmA; // (T+1)*(J+1)
  uint8_t *opB, *pmB;
  size_t T, J;
} Backptr;
static inline size_t IDX(size_t t,size_t j,size_t Jp1){ return t*Jp1 + j; }
static Backptr backptr_alloc(size_t T,size_t J){
  Backptr bp; bp.T=T; bp.J=J;
  size_t cells=(T+1)*(J+1);
  bp.opA=xcalloc(cells,1); bp.pmA=xcalloc(cells,1);
  bp.opB=xcalloc(cells,1); bp.pmB=xcalloc(cells,1);
  return bp;
}
static void backptr_free(Backptr *bp){
  free(bp->opA); free(bp->pmA); free(bp->opB); free(bp->pmB);
  bp->opA=bp->pmA=bp->opB=bp->pmB=NULL;
}

// ---------- Precompute helpers ----------

// gap_any[t]=1 if A[t-1]=='-' or B[t-1]=='-'
static void compute_gap_any(const char *Aaln,const char *Baln,size_t T,uint8_t *gap_any){
  for(size_t t=1;t<=T;t++){ gap_any[t]=(Aaln[t-1]=='-' || Baln[t-1]=='-')?1:0; }
}
static void prefix_u8(const uint8_t *arr,size_t n,int *pref){ // 1..n
  pref[0]=0;
  for(size_t i=1;i<=n;i++) pref[i]=pref[i-1]+(int)arr[i];
}
static inline int any_in_window_pref(const int *pref,size_t n,int l,int r){
  if(l<1) l=1; if(r>(int)n) r=(int)n;
  if(l>r) return 0;
  return (pref[r]-pref[l-1])>0;
}

// informative on aligned t: A!=B and both non-gap
static void compute_inf_aligned(const char *Aaln,const char *Baln,size_t T,uint8_t *inf_t){
  for(size_t t=1;t<=T;t++){
    char a=Aaln[t-1], b=Baln[t-1];
    inf_t[t] = (a!='-' && b!='-' && a!=b) ? 1 : 0;
  }
}

// informative on ungapped index u: A*!=B*
static void compute_inf_ungapped(const char *Astar,const char *Bstar,size_t Lmin,uint8_t *inf_u){
  for(size_t u=1;u<=Lmin;u++){
    char a=Astar[u-1], b=Bstar[u-1];
    inf_u[u] = (a!=b)?1:0;
  }
}

// unary bias betaA/B[j] from windowed informative matches on ungapped coordinates
static void compute_beta(const char *Astar,const char *Bstar,size_t Lmin,
                         const char *Y,size_t J,int W,int kappa,
                         const uint8_t *inf_u,int *betaA,int *betaB)
{
  int half=W/2;
  size_t Nmin = (J<Lmin? J : Lmin);

  // wA[u]=1 if Y[u]==A*[u] and inf[u]==1 else 0, similarly wB
  int *prefA = xcalloc(Nmin+1, sizeof(int));
  int *prefB = xcalloc(Nmin+1, sizeof(int));
  for(size_t u=1; u<=Nmin; ++u){
    int wA = (Y[u-1]==Astar[u-1] && inf_u[u]) ? 1 : 0;
    int wB = (Y[u-1]==Bstar[u-1] && inf_u[u]) ? 1 : 0;
    prefA[u] = prefA[u-1] + wA;
    prefB[u] = prefB[u-1] + wB;
  }

  for(size_t j=1;j<=J;j++){
    if(j<=Nmin){
      int L=(int)j - half, R=(int)j + half;
      if(L<1) L=1; if(R>(int)Nmin) R=(int)Nmin;
      int SA = (L<=R)? (prefA[R]-prefA[L-1]) : 0;
      int SB = (L<=R)? (prefB[R]-prefB[L-1]) : 0;
      int s = SA - SB;
      betaA[j] = -kappa * s;
      betaB[j] = +kappa * s;
    }else{
      betaA[j]=0; betaB[j]=0;
    }
  }
  free(prefA); free(prefB);
}

// ---------- Core DP ----------
typedef struct {
  int total, mismatches, insertions, deletions, switches;
} Counts;

static Counts red_align(
  const char *Aaln,const char *Baln,size_t T,
  const char *Astar,const char *Bstar,size_t Lmin,
  const char *Y,size_t J,
  int mu_inf, int mu_eq, int reward,
  int gap, int lambda_open, int stick,
  int gap_swell, int gap_near,
  const uint8_t *inf_t_aligned,
  const int *betaA, const int *betaB,
  const int *prefGap // size T (prefix over gap_any) to query nearGap
){
  const size_t Jp1 = J+1;
  int *prevA = xmalloc(Jp1*sizeof(int));
  int *prevB = xmalloc(Jp1*sizeof(int));
  int *curA  = xmalloc(Jp1*sizeof(int));
  int *curB  = xmalloc(Jp1*sizeof(int));

  Backptr bp = backptr_alloc(T,J);

  // init (0,0)
  prevA[0]=0; prevB[0]=0;

  // row 0: insertions include beta at each j
  for(size_t j=1;j<=J;j++){
    prevA[j] = prevA[j-1] + gap + betaA[j];
    prevB[j] = prevB[j-1] + gap + betaB[j];
    size_t idx=IDX(0,j,Jp1);
    bp.opA[idx]=OP_INS; bp.pmA[idx]=0;
    bp.opB[idx]=OP_INS; bp.pmB[idx]=1;
  }

  // DP over t rows
  for(size_t t=1;t<=T;t++){
    // j=0: deletions (no beta)
    curA[0] = prevA[0] + gap; size_t idxA0=IDX(t,0,Jp1); bp.opA[idxA0]=OP_DEL; bp.pmA[idxA0]=0;
    curB[0] = prevB[0] + gap; size_t idxB0=IDX(t,0,Jp1); bp.opB[idxB0]=OP_DEL; bp.pmB[idxB0]=1;

    char a=Aaln[t-1], b=Baln[t-1];
    int mu_t = (inf_t_aligned[t] ? mu_inf : mu_eq);

    // near-gap switch swell
    int lam_t = lambda_open;
    // gap_any prefix: any gap in [t-gap_near, t+gap_near] ?
    if(gap_near>0 && prefGap){
      int L=(int)t - gap_near, R=(int)t + gap_near;
      if(L<1) L=1; if(R>(int)T) R=(int)T;
      if((prefGap[R]-prefGap[L-1])>0) lam_t += gap_swell;
    }

    for(size_t j=1;j<=J;j++){
      char y = Y[j-1];

      // Mode A
      {
        int best=INT_MAX, op=OP_NONE, pm=0;

        if(a!='-'){
          int sub = (a==y)? 0 : mu_t;
          if(inf_t_aligned[t] && a==y && reward>0) sub -= reward;
          int cand = prevA[j-1] + sub - stick + betaA[j]; // diag stay
          if(cand<best){ best=cand; op=OP_DIAG; pm=0; }

          cand = prevB[j-1] + lam_t + sub + betaA[j]; // switch
          if(cand<best){ best=cand; op=OP_SWITCH; pm=1; }
        }
        // insertion (consume child)
        int cand = curA[j-1] + gap + betaA[j];
        if(cand<best){ best=cand; op=OP_INS; pm=0; }

        // deletion (consume parent)
        cand = prevA[j] + gap;
        if(cand<best){ best=cand; op=OP_DEL; pm=0; }

        curA[j]=best;
        size_t idx=IDX(t,j,Jp1);
        bp.opA[idx]=(uint8_t)op; bp.pmA[idx]=(uint8_t)pm;
      }

      // Mode B
      {
        int best=INT_MAX, op=OP_NONE, pm=1;

        if(b!='-'){
          int sub = (b==y)? 0 : mu_t;
          if(inf_t_aligned[t] && b==y && reward>0) sub -= reward;
          int cand = prevB[j-1] + sub - stick + betaB[j]; // diag stay
          if(cand<best){ best=cand; op=OP_DIAG; pm=1; }

          cand = prevA[j-1] + lam_t + sub + betaB[j]; // switch
          if(cand<best){ best=cand; op=OP_SWITCH; pm=0; }
        }
        // insertion
        int cand = curB[j-1] + gap + betaB[j];
        if(cand<best){ best=cand; op=OP_INS; pm=1; }

        // deletion
        cand = prevB[j] + gap;
        if(cand<best){ best=cand; op=OP_DEL; pm=1; }

        curB[j]=best;
        size_t idx=IDX(t,j,Jp1);
        bp.opB[idx]=(uint8_t)op; bp.pmB[idx]=(uint8_t)pm;
      }
    }

    // rotate
    int *tmpA=prevA; prevA=curA; curA=tmpA;
    int *tmpB=prevB; prevB=curB; curB=tmpB;
  }

  // terminal
  int endA=prevA[J], endB=prevB[J];
  int mstar = (endA<=endB)? 0:1;
  int total = (mstar==0)? endA : endB;

  // backtrace
  Counts cnt = { .total=total, .mismatches=0, .insertions=0, .deletions=0, .switches=0 };
  char *orig = NULL; // fill later (we only count during backtrace here)
  // For output, we need origins per child; but we return counts here.
  // We'll backtrace in the caller where we have 'orig' buffer.
  // (To keep a single function interface, we do backtrace here with a dummy pass.)
  // Instead, we expose bp and arrays? Simpler: do the real backtrace now:
  // (Allocate a temp to store origins, then copy back outside via global? Better: refactor.)
  // Easier: We'll do a mini backtrace here to count, then caller redoes backtrace to fill origins.
  // But double backtrace is fine (O(T+J)). Let's just return counts after a true backtrace performed below in caller.

  // cleanup temps; backptr stays needed; but we'll return bp etc.
  // To keep API simple, we'll duplicate minimal backtrace code outside using bp.
  // Free here and recompute would be wasteful. So change plan: we will perform full backtrace here,
  // but we need 'orig' from caller. Let's change function signature to accept 'orig' and fill it.

  // -- unreachable --
  free(prevA); free(prevB); free(curA); free(curB);
  backptr_free((Backptr*)&bp);
  return cnt;
}

// We'll implement a wrapper that runs DP and backtrace, filling origins and counts.
typedef struct {
  int total, mismatches, insertions, deletions, switches;
  int endmode; // 0=A,1=B
} Result;

static Result red_run(
  const char *Aaln,const char *Baln,size_t T,
  const char *Astar,const char *Bstar,size_t Lmin,
  const char *Y,size_t J,
  int mu_inf, int mu_eq, int reward,
  int gap, int lambda_open, int stick,
  int gap_swell, int gap_near,
  const uint8_t *inf_t_aligned,
  const int *betaA, const int *betaB,
  const int *prefGap,
  char *orig // output labels length J
){
  const size_t Jp1 = J+1;
  int *prevA = xmalloc((Jp1)*sizeof(int));
  int *prevB = xmalloc((Jp1)*sizeof(int));
  int *curA  = xmalloc((Jp1)*sizeof(int));
  int *curB  = xmalloc((Jp1)*sizeof(int));

  Backptr bp = backptr_alloc(T,J);

  prevA[0]=0; prevB[0]=0;
  for(size_t j=1;j<=J;j++){
    prevA[j] = prevA[j-1] + gap + betaA[j];
    prevB[j] = prevB[j-1] + gap + betaB[j];
    size_t idx=IDX(0,j,Jp1);
    bp.opA[idx]=OP_INS; bp.pmA[idx]=0;
    bp.opB[idx]=OP_INS; bp.pmB[idx]=1;
  }

  uint8_t *inf_t = (uint8_t*)inf_t_aligned; // alias
  for(size_t t=1;t<=T;t++){
    curA[0] = prevA[0] + gap; size_t idxA0=IDX(t,0,Jp1); bp.opA[idxA0]=OP_DEL; bp.pmA[idxA0]=0;
    curB[0] = prevB[0] + gap; size_t idxB0=IDX(t,0,Jp1); bp.opB[idxB0]=OP_DEL; bp.pmB[idxB0]=1;

    char a=Aaln[t-1], b=Baln[t-1];
    int mu_t = (inf_t[t] ? mu_inf : mu_eq);
    int lam_t = lambda_open;
    if(gap_near>0 && prefGap){
      int L=(int)t-gap_near, R=(int)t+gap_near;
      if(L<1) L=1; if(R>(int)T) R=(int)T;
      if((prefGap[R]-prefGap[L-1])>0) lam_t += gap_swell;
    }

    for(size_t j=1;j<=J;j++){
      char y=Y[j-1];

      // A mode
      {
        int best=INT_MAX, op=OP_NONE, pm=0;
        if(a!='-'){
          int sub=(a==y)?0:mu_t;
          if(inf_t[t] && a==y && reward>0) sub -= reward;
          int cand=prevA[j-1]+sub - stick + betaA[j];
          if(cand<best){best=cand;op=OP_DIAG;pm=0;}
          cand=prevB[j-1]+lam_t+sub + betaA[j];
          if(cand<best){best=cand;op=OP_SWITCH;pm=1;}
        }
        int cand=curA[j-1]+gap + betaA[j];
        if(cand<best){best=cand;op=OP_INS;pm=0;}
        cand=prevA[j]+gap;
        if(cand<best){best=cand;op=OP_DEL;pm=0;}
        curA[j]=best;
        size_t idx=IDX(t,j,Jp1); bp.opA[idx]=(uint8_t)op; bp.pmA[idx]=(uint8_t)pm;
      }

      // B mode
      {
        int best=INT_MAX, op=OP_NONE, pm=1;
        if(b!='-'){
          int sub=(b==y)?0:mu_t;
          if(inf_t[t] && b==y && reward>0) sub -= reward;
          int cand=prevB[j-1]+sub - stick + betaB[j];
          if(cand<best){best=cand;op=OP_DIAG;pm=1;}
          cand=prevA[j-1]+lam_t+sub + betaB[j];
          if(cand<best){best=cand;op=OP_SWITCH;pm=0;}
        }
        int cand=curB[j-1]+gap + betaB[j];
        if(cand<best){best=cand;op=OP_INS;pm=1;}
        cand=prevB[j]+gap;
        if(cand<best){best=cand;op=OP_DEL;pm=1;}
        curB[j]=best;
        size_t idx=IDX(t,j,Jp1); bp.opB[idx]=(uint8_t)op; bp.pmB[idx]=(uint8_t)pm;
      }
    }

    int *tmpA=prevA; prevA=curA; curA=tmpA;
    int *tmpB=prevB; prevB=curB; curB=tmpB;
  }

  int endA=prevA[J], endB=prevB[J];
  int endmode = (endA<=endB)? 0:1;
  int total = (endmode==0)? endA : endB;

  // backtrace, fill orig + counts
  Result R={0,0,0,0,0,endmode}; R.total=total;
  size_t t=T, j=J; int m=endmode;
  while(t>0 || j>0){
    size_t idx=IDX(t,j,Jp1);
    int op, pm; char pchar=0;
    if(m==0){ op=bp.opA[idx]; pm=bp.pmA[idx]; if(t>0) pchar=Aaln[t-1]; }
    else    { op=bp.opB[idx]; pm=bp.pmB[idx]; if(t>0) pchar=Baln[t-1]; }

    if(op==OP_DIAG){
      if(j>0 && t>0){
        if(pchar!=Y[j-1]) R.mismatches++;
        orig[j-1] = (m==0)? 'A':'B';
        t--; j--; m=pm;
      } else break;
    } else if(op==OP_SWITCH){
      if(j>0 && t>0){
        if(pchar!=Y[j-1]) R.mismatches++;
        orig[j-1] = (m==0)? 'A':'B';
        R.switches++;
        t--; j--; m=pm;
      } else break;
    } else if(op==OP_INS){
      if(j>0){
        orig[j-1]='I';
        R.insertions++;
        j--;
      } else break;
    } else if(op==OP_DEL){
      if(t>0){
        R.deletions++;
        t--;
      } else break;
    } else {
      if(j>0){ orig[j-1]='I'; j--; }
      else if(t>0){ t--; }
      else break;
    }
  }
  // fill remaining undefined
  char last='A';
  for(size_t k=0;k<J;k++){
    if(orig[k]!='A' && orig[k]!='B' && orig[k]!='I') orig[k]=last;
    if(orig[k]!='I') last=orig[k];
  }

  backptr_free(&bp);
  free(prevA); free(prevB); free(curA); free(curB);
  return R;
}

// ---------- main ----------
int main(int argc, char **argv){
  const char *parents_path=get_arg(argc,argv,"--parents-aln");
  const char *child_path  =get_arg(argc,argv,"--child");
  const char *outprefix   =get_arg(argc,argv,"--out");
  if(!parents_path || !child_path || !outprefix){ usage(); return 2; }

  // legacy params
  int mismatch = get_arg(argc,argv,"--mismatch") ? atoi(get_arg(argc,argv,"--mismatch")) : 1;
  int gap      = get_arg(argc,argv,"--gap")      ? atoi(get_arg(argc,argv,"--gap"))      : 2;
  int lambda_open = get_arg(argc,argv,"--switch")? atoi(get_arg(argc,argv,"--switch"))    : 1;

  // RED params
  int W        = get_arg(argc,argv,"--win")      ? atoi(get_arg(argc,argv,"--win"))      : 21;
  int kappa    = get_arg(argc,argv,"--kappa")    ? atoi(get_arg(argc,argv,"--kappa"))    : 1;
  int stick    = get_arg(argc,argv,"--stick")    ? atoi(get_arg(argc,argv,"--stick"))    : 0;
  int mu_inf   = get_arg(argc,argv,"--mu-inf")   ? atoi(get_arg(argc,argv,"--mu-inf"))   : -1;
  int mu_eq    = get_arg(argc,argv,"--mu-eq")    ? atoi(get_arg(argc,argv,"--mu-eq"))    : -1;
  int gap_swell= get_arg(argc,argv,"--gap-swell")? atoi(get_arg(argc,argv,"--gap-swell")): 1;
  int gap_near = get_arg(argc,argv,"--gap-near") ? atoi(get_arg(argc,argv,"--gap-near")) : 10;
  int reward   = get_arg(argc,argv,"--reward")   ? atoi(get_arg(argc,argv,"--reward"))   : 0;

  if(W<=1) W=21; if((W%2)==0) W+=1;
  if(kappa<0) kappa=0;
  if(stick<0) stick=0;
  if(reward<0) reward=0;

  // defaults for μ
  if(mu_eq<0)  mu_eq  = mismatch;
  if(mu_inf<0) mu_inf = (mismatch>0)? (2*mismatch) : 0; // stronger penalty on informative sites

  // Read FASTAs
  FastaVec P; read_fasta(parents_path,&P);
  if(P.n!=2) die("parents-aln must contain exactly two sequences");
  if(P.recs[0].len != P.recs[1].len) die("parents-aln sequences must have equal length");
  const char *Aaln=P.recs[0].seq, *Baln=P.recs[1].seq; size_t T=P.recs[0].len;

  FastaVec C; read_fasta(child_path,&C);
  if(C.n==0) die("no child sequences found");

  // Ungap parents for bootscan bias
  size_t LA=0, LB=0;
  char *Astar=ungap(Aaln,&LA), *Bstar=ungap(Baln,&LB);
  size_t Lmin = (LA<LB?LA:LB);

  // Precompute aligned informative and near-gap prefix
  uint8_t *inf_t = xcalloc(T+1,1);
  compute_inf_aligned(Aaln,Baln,T,inf_t);

  uint8_t *gap_any = xcalloc(T+1,1);
  compute_gap_any(Aaln,Baln,T,gap_any);
  int *prefGap = xcalloc(T+1, sizeof(int));
  prefix_u8(gap_any, T, prefGap);

  // For each child
  for(size_t ci=0; ci<C.n; ++ci){
    const char *Y=C.recs[ci].seq; size_t J=C.recs[ci].len;

    // Compute ungapped-informative and betaA/B[j]
    size_t Nmin = (J<Lmin? J : Lmin);
    uint8_t *inf_u = xcalloc((Nmin>0?Nmin:1)+1,1);
    compute_inf_ungapped(Astar,Bstar,Nmin,inf_u);

    int *betaA = xcalloc(J+1,sizeof(int));
    int *betaB = xcalloc(J+1,sizeof(int));
    compute_beta(Astar,Bstar,Lmin, Y,J, W,kappa, inf_u, betaA,betaB);

    char *orig = xmalloc(J);
    Result R = red_run(
      Aaln,Baln,T, Astar,Bstar,Lmin, Y,J,
      mu_inf,mu_eq,reward, gap, lambda_open, stick,
      gap_swell, gap_near, inf_t, betaA,betaB, prefGap, orig
    );

    // write outputs
    char path1[4096], path2[4096];
    // C.recs[ci].name -> ""
    snprintf(path1,sizeof(path1), "%s.%s.origins.tsv",  outprefix, "red");
    snprintf(path2,sizeof(path2), "%s.%s.summary.tsv",  outprefix, "red");
    write_origins_tsv(path1, orig, J);
    write_summary_tsv(path2, R.total, R.mismatches, R.insertions, R.deletions, R.switches);

    fprintf(stderr,"[OK] red %s  total=%d  mis=%d ins=%d del=%d sw=%d  (end=%c)  W=%d kappa=%d stick=%d mu_inf=%d mu_eq=%d\n",
            C.recs[ci].name, R.total, R.mismatches, R.insertions, R.deletions, R.switches,
            (R.endmode? 'B':'A'), W,kappa,stick,mu_inf,mu_eq);

    free(orig); free(inf_u); free(betaA); free(betaB);
  }

  // cleanup
  for(size_t i=0;i<P.n;i++){ free(P.recs[i].name); free(P.recs[i].seq); }
  for(size_t i=0;i<C.n;i++){ free(C.recs[i].name); free(C.recs[i].seq); }
  free(P.recs); free(C.recs);
  free(Astar); free(Bstar);
  free(inf_t); free(gap_any); free(prefGap);
  return 0;
}

