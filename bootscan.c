// bootscan.c — Sliding-window bootscan (C, single-file CLI)
// Author: Katsuhiko Hayashi
// Build:   gcc -O3 -std=c11 -o bootscan bootscan.c
// Usage:   ./bootscan --parents-aln parents.aligned.fasta --child child.fasta \
//                     --out out/prefix [--window 101] [--step 1]
//
// Output per child:
//   <outprefix>.<childName>.origins.tsv   (j \t parent[A|B|INS]; )
//   <outprefix>.<childName>.summary.tsv   (total,mismatches,insertions=0,deletions=0,switches)
//
// Notes:
// - Parents: exactly two aligned sequences (equal length), may contain '-'.
// - Child:   one or more ungapped sequences.
// - This is a lightweight alternative to RED: it assigns A/B per position by
//   comparing windowed Hamming mismatches to ungapped parents A* and B*.

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

typedef struct { char *name; char *seq; size_t len; } FastaRec;
typedef struct { FastaRec *recs; size_t n, cap; } FastaVec;

static void die(const char *msg){ fprintf(stderr, "[ERROR] %s\n", msg); exit(1); }
static void *xmalloc(size_t n){ void *p=malloc(n); if(!p) die("malloc"); return p; }
static void *xrealloc(void *q, size_t n){ void *p=realloc(q,n); if(!p) die("realloc"); return p; }
static char *xstrdup(const char *s){ char *p=strdup(s); if(!p) die("strdup"); return p; }

static void to_upper_inplace(char *s){ for(;*s;++s)*s=(char)toupper((unsigned char)*s); }

static void fasta_vec_init(FastaVec *v){ v->recs=NULL; v->n=0; v->cap=0; }
static void fasta_vec_push(FastaVec *v, FastaRec rec){
  if(v->n==v->cap){ v->cap = v->cap? v->cap*2 : 8; v->recs = xrealloc(v->recs, v->cap*sizeof(FastaRec)); }
  v->recs[v->n++] = rec;
}

// Read multi-FASTA (concatenate wrapped lines; strip whitespace)
static void read_fasta(const char *path, FastaVec *out){
  fasta_vec_init(out);
  FILE *fp = fopen(path,"r"); if(!fp){ perror(path); die("open fasta"); }
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
      if(!buf){ bcap = (size_t)n+64; buf = xmalloc(bcap); buf[0]=0; }
      if(blen+(size_t)n+1 > bcap){ bcap=(blen+(size_t)n+1)*2; buf=xrealloc(buf,bcap); }
      for(size_t i=0;i<(size_t)n;i++){ if(!isspace((unsigned char)line[i])) buf[blen++]=line[i]; }
      buf[blen]=0;
    }
  }
  if(name){
    FastaRec rec={ xstrdup(name), xstrdup(buf?buf:""), blen };
    fasta_vec_push(out, rec);
  }
  free(line); if(buf) free(buf); fclose(fp);
  for(size_t i=0;i<out->n;i++){ to_upper_inplace(out->recs[i].seq); out->recs[i].len=strlen(out->recs[i].seq); }
}

static char *ungap(const char *s, size_t *out_len){
  size_t n=strlen(s);
  char *u=xmalloc(n+1); size_t k=0;
  for(size_t i=0;i<n;i++){ if(s[i]!='-') u[k++]=s[i]; }
  u[k]=0; if(out_len) *out_len=k; return u;
}

static const char *get_arg(int argc, char **argv, const char *key){
  for(int i=1;i<argc;i++){ if(strcmp(argv[i],key)==0){ if(i+1<argc && argv[i+1][0]!='-') return argv[i+1]; else return "true"; } }
  return NULL;
}

static void usage(){
  fprintf(stderr,
    "Usage: bootscan --parents-aln PARENTS.fasta --child CHILD.fasta --out OUTPREFIX "
    "[--window 101] [--step 1]\n");
}

static void write_origins_tsv(const char *path, const char *orig, size_t J){
  FILE *fp=fopen(path,"w"); if(!fp){ perror(path); die("open origins out"); }
  fprintf(fp,"j\tparent\n");
  for(size_t j=1;j<=J;j++){
    char c=orig[j-1];
    if(c=='I') fprintf(fp,"%zu\tINS\n", j);           // not used here, but kept for compatibility
    else       fprintf(fp,"%zu\t%c\n",  j, c);        // 'A' or 'B'
  }
  fclose(fp);
}

static void write_summary_tsv(const char *path, int total, int mism, int switches){
  FILE *fp=fopen(path,"w"); if(!fp){ perror(path); die("open summary out"); }
  fprintf(fp,"total\tmismatches\tinsertions\tdeletions\tswitches\n");
  fprintf(fp,"%d\t%d\t%d\t%d\t%d\n", total, mism, 0, 0, switches);
  fclose(fp);
}

// count mismatches on [L..R] (1-based inclusive, clipped to [1..Nmin])
static inline int window_mismatches(const char *Y, const char *P, size_t Nmin, long L, long R){
  if(L<1) L=1; if(R>(long)Nmin) R=(long)Nmin;
  int mm=0;
  for(long i=L;i<=R;i++){ char y=Y[i-1], p=P[i-1]; if(y!=p) mm++; }
  return mm;
}

int main(int argc, char **argv){
  const char *parents_path = get_arg(argc, argv, "--parents-aln");
  const char *child_path   = get_arg(argc, argv, "--child");
  const char *outprefix    = get_arg(argc, argv, "--out");
  const char *win_s        = get_arg(argc, argv, "--window");
  const char *step_s       = get_arg(argc, argv, "--step");

  if(!parents_path || !child_path || !outprefix){ usage(); return 2; }
  int W = win_s ? atoi(win_s) : 101; if(W<=1) W=101; if((W%2)==0) W+=1; // make odd
  int STEP = step_s ? atoi(step_s) : 1; if(STEP<1) STEP=1;

  // read parents (aligned, length T with '-')
  FastaVec P; read_fasta(parents_path, &P);
  if(P.n!=2) die("parents-aln must contain exactly two sequences");
  if(P.recs[0].len != P.recs[1].len) die("parents-aln sequences must have equal length");
  // ungap parents -> A*, B*
  size_t LA=0, LB=0;
  char *Astar = ungap(P.recs[0].seq, &LA);
  char *Bstar = ungap(P.recs[1].seq, &LB);
  size_t Lmin = (LA<LB?LA:LB);

  // read child
  FastaVec C; read_fasta(child_path, &C);
  if(C.n==0) die("no child sequences found");

  for(size_t ci=0; ci<C.n; ++ci){
    const char *Y = C.recs[ci].seq; size_t J = C.recs[ci].len;
    size_t Nmin = (J<Lmin? J : Lmin);
    if(J!=Lmin){
      fprintf(stderr,"[WARN] length mismatch child=%zu vs parents(min)=%zu; clipping to %zu\n", J, Lmin, Nmin);
    }

    char *orig = xmalloc(J);
    // default previous label 'A' for ties / leading positions
    char prev = 'A';
    int total_mism=0, switches=0;

    int half = W/2;
    // assign label at every position j (1..J), comparing only up to Nmin
    for(size_t j=1;j<=J;j++){
      char label = prev;
      if(j<=Nmin){
        long L = (long)j - half;
        long R = (long)j + half;
        int mmA = window_mismatches(Y, Astar, Nmin, L, R);
        int mmB = window_mismatches(Y, Bstar, Nmin, L, R);
        if(mmA < mmB) label='A';
        else if(mmB < mmA) label='B';
        else label=prev; // tie -> keep previous
        // per-base mismatch vs chosen parent (for summary)
        if(label=='A'){ if(Y[j-1]!=Astar[j-1]) total_mism++; }
        else          { if(Y[j-1]!=Bstar[j-1]) total_mism++; }
      }else{
        // beyond Nmin: keep previous label
        label = prev;
      }
      orig[j-1] = label;
      if(j>1 && label!=prev) switches++;
      prev = label;

      // STEP>1 の場合でも origins は全位置に出す（スムージング強度は窓で表現）
      // もし本当に subsampling したいなら、ここで j+=STEP-1 する。ただし評価互換性のため全点割当てにしています。
    }

    // total: use total mismatches as a simple cost proxy
    int total_cost = total_mism;

    // write outputs
    char path1[4096], path2[4096];
    // C.recs[ci].name -> ""
    snprintf(path1,sizeof(path1), "%s.%s.origins.tsv",  outprefix, "bs");
    snprintf(path2,sizeof(path2), "%s.%s.summary.tsv",  outprefix, "bs");
    write_origins_tsv(path1, orig, J);
    write_summary_tsv(path2, total_cost, total_mism, switches);

    fprintf(stderr,"[OK] bootscan %s  W=%d  switches=%d  mism=%d\n",
            C.recs[ci].name, W, switches, total_mism);

    free(orig);
  }

  // cleanup
  for(size_t i=0;i<P.n;i++){ free(P.recs[i].name); free(P.recs[i].seq); }
  for(size_t i=0;i<C.n;i++){ free(C.recs[i].name); free(C.recs[i].seq); }
  free(P.recs); free(C.recs);
  free(Astar); free(Bstar);
  return 0;
}

