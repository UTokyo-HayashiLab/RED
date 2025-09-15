// red_band.c : RED-B with Ukkonen-style banding (unary-bias guided)
// Author: Katsuhiko Hayashi
// CLI: same flags as red.c (+ --band). Supports --mismatch (sets mu-eq = mu-inf).
// Outputs:
//   res.<child>.origins.tsv  (headerless: "j\tA|B")
//   res.<child>.summary.tsv  (2-line TSV):
//     total\tmismatches\tinsertions\tdeletions\tswitches
//     <int>\t<int>\t<int>\t<int>\t<int>
//
// Build: gcc -O3 -std=c11 -o red_band red_band.c

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>

#define MAX_NAME 256
#define MAX_LINE 8192

static void die(const char*fmt, ...) {
  va_list ap; va_start(ap, fmt);
  vfprintf(stderr, fmt, ap); fputc('\n', stderr);
  va_end(ap);
  exit(1);
}
static void to_upper_str(char*s){ for(;*s;++s) *s=(char)toupper((unsigned char)*s); }
static inline int is_base(char c){ return (c=='A'||c=='C'||c=='G'||c=='T'); }
static inline int imax(int a,int b){ return a>b?a:b; }
static inline int imin(int a,int b){ return a<b?a:b; }
static inline int iround(double x){ return (int)llround(x); }

typedef struct {
  int mismatch_eq;   // mu_eq
  int mismatch_inf;  // mu_inf
  int reward;        // informative match reward (>=0)
  int gap;           // insertion/deletion cost
  int sw_base;       // lambda (switch cost)
  int gap_swell;     // near-gap extra for switch
  int gap_near;      // +/-g for near-gap
  int stick;         // stay stickiness reward
  int win;           // unary window (odd)
  int kappa;         // unary scaling
  int band;          // band half-width
} Params;

static void set_default_params(Params *p){
  p->mismatch_eq = 1;
  p->mismatch_inf= 2;
  p->reward      = 0;
  p->gap         = 1;
  p->sw_base     = 1;
  p->gap_swell   = 1;
  p->gap_near    = 10;
  p->stick       = 1;
  p->win         = 21;
  p->kappa       = 1;
  p->band        = -1; // auto: max(10, W)
}

typedef struct {
  char name[MAX_NAME];
  char *seq; // 1-based
  int len;
} Fasta;

static Fasta read_fasta_one(const char*path){
  FILE*fp=fopen(path,"r"); if(!fp) die("open %s: %s", path, strerror(errno));
  char line[MAX_LINE];
  Fasta f; memset(&f,0,sizeof(f));
  char *buf=NULL; size_t bufsz=0, used=0; bool seen=false;

  while(fgets(line,sizeof(line),fp)){
    if(line[0]=='>'){
      if(seen) break;
      seen=true;
      char *p=line+1; while(*p && isspace((unsigned char)*p)) ++p;
      char *q=p; while(*q && !isspace((unsigned char)*q)) ++q;
      size_t n=(size_t)(q-p); if(n>=MAX_NAME) n=MAX_NAME-1;
      memcpy(f.name,p,n); f.name[n]='\0';
    }else{
      for(char *p=line; *p; ++p){
        if(isspace((unsigned char)*p)) continue;
        if(bufsz-used<2){ bufsz = bufsz? bufsz*2 : 1024; buf=realloc(buf,bufsz); if(!buf) die("OOM"); }
        buf[used++] = *p;
      }
    }
  }
  fclose(fp);
  if(!seen) die("FASTA %s: no header", path);
  buf[used]='\0'; to_upper_str(buf);
  f.len=(int)used;
  f.seq=(char*)malloc((size_t)f.len+2); if(!f.seq) die("OOM");
  f.seq[0]='?'; memcpy(f.seq+1, buf, (size_t)used+1);
  free(buf);
  return f;
}

typedef struct { Fasta A,B; } Parents;

static Parents read_parents_aln(const char*path){
  FILE*fp=fopen(path,"r"); if(!fp) die("open %s: %s", path, strerror(errno));
  char line[MAX_LINE];
  Fasta recs[2]; int r=0;
  char *buf=NULL; size_t bufsz=0, used=0;

  while(fgets(line,sizeof(line),fp)){
    if(line[0]=='>'){
      if(r==2) break;
      if(r>0){
        buf[used]='\0'; to_upper_str(buf);
        recs[r-1].len=(int)used;
        recs[r-1].seq=(char*)malloc((size_t)used+2);
        if(!recs[r-1].seq) die("OOM");
        recs[r-1].seq[0]='?';
        memcpy(recs[r-1].seq+1, buf, (size_t)used+1);
        used=0;
      }
      ++r;
      char *p=line+1; while(*p && isspace((unsigned char)*p)) ++p;
      char *q=p; while(*q && !isspace((unsigned char)*q)) ++q;
      size_t n=(size_t)(q-p); if(n>=MAX_NAME) n=MAX_NAME-1;
      memcpy(recs[r-1].name,p,n); recs[r-1].name[n]='\0';
    }else{
      for(char *p=line; *p; ++p){
        if(isspace((unsigned char)*p)) continue;
        if(bufsz-used<2){ bufsz=bufsz?bufsz*2:2048; buf=realloc(buf,bufsz); if(!buf) die("OOM"); }
        buf[used++]=*p;
      }
    }
  }
  if(r<2) die("parents FASTA must have exactly 2 records");
  buf[used]='\0'; to_upper_str(buf);
  recs[r-1].len=(int)used;
  recs[r-1].seq=(char*)malloc((size_t)used+2);
  if(!recs[r-1].seq) die("OOM");
  recs[r-1].seq[0]='?';
  memcpy(recs[r-1].seq+1, buf, (size_t)used+1);
  free(buf);
  fclose(fp);
  if(recs[0].len!=recs[1].len) die("parent alignment length mismatch: %d vs %d", recs[0].len, recs[1].len);
  Parents P; P.A=recs[0]; P.B=recs[1]; return P;
}

typedef struct {
  int T, J;
  char *A,*B; // aligned parents (1..T)
  char *Y;    // child (1..J)

  uint8_t *I;
  int *mu;
  uint8_t *nearGap;

  int *betaA, *betaB;

  int *mapA, *mapB; int La, Lb;

  int *tHatA, *tHatB;
  int *loA, *hiA, *loB, *hiB;
} Pre;

static void build_informative_and_neargap(Pre *pr, const Params *pm){
  pr->I = (uint8_t*)calloc((size_t)pr->T+1,1);
  pr->mu= (int*)     calloc((size_t)pr->T+1,sizeof(int));
  pr->nearGap=(uint8_t*)calloc((size_t)pr->T+1,1);
  if(!pr->I||!pr->mu||!pr->nearGap) die("OOM");

  for(int t=1;t<=pr->T;++t){
    char a=pr->A[t], b=pr->B[t];
    int inf = (is_base(a)&&is_base(b)&&(a!=b)) ? 1:0;
    pr->I[t]=(uint8_t)inf;
    pr->mu[t]= inf? pm->mismatch_inf : pm->mismatch_eq;
  }
  int g=pm->gap_near;
  for(int t=1;t<=pr->T;++t){
    int near=0;
    int L = (t-g>1)? t-g:1;
    int R = (t+g<pr->T)? t+g:pr->T;
    for(int q=L;q<=R;++q){
      if(pr->A[q]=='-' || pr->B[q]=='-'){ near=1; break; }
    }
    pr->nearGap[t]=(uint8_t)near;
  }
}

static void build_ungap_maps(Pre*pr){
  pr->mapA=(int*)malloc(sizeof(int)*(pr->T+1));
  pr->mapB=(int*)malloc(sizeof(int)*(pr->T+1));
  if(!pr->mapA||!pr->mapB) die("OOM");
  pr->La=0; pr->Lb=0;
  for(int t=1;t<=pr->T;++t) if(is_base(pr->A[t])) pr->mapA[++pr->La]=t;
  for(int t=1;t<=pr->T;++t) if(is_base(pr->B[t])) pr->mapB[++pr->Lb]=t;
}

static void build_unary_bias(Pre*pr, const Params*pm){
  char *Au=(char*)malloc((size_t)pr->T+1);
  char *Bu=(char*)malloc((size_t)pr->T+1);
  if(!Au||!Bu) die("OOM");
  int La=0, Lb=0;
  for(int t=1;t<=pr->T;++t) if(is_base(pr->A[t])) Au[++La]=pr->A[t];
  for(int t=1;t<=pr->T;++t) if(is_base(pr->B[t])) Bu[++Lb]=pr->B[t];
  int L = imin(imin(La,Lb), pr->J);

  int *MA=(int*)calloc((size_t)pr->J+1,sizeof(int));
  int *MB=(int*)calloc((size_t)pr->J+1,sizeof(int));
  if(!MA||!MB) die("OOM");
  for(int j=1;j<=pr->J;++j){
    int mA=(j<=L && pr->Y[j]==Au[j])?1:0;
    int mB=(j<=L && pr->Y[j]==Bu[j])?1:0;
    MA[j]=MA[j-1]+mA;
    MB[j]=MB[j-1]+mB;
  }

  pr->betaA=(int*)calloc((size_t)pr->J+1,sizeof(int));
  pr->betaB=(int*)calloc((size_t)pr->J+1,sizeof(int));
  if(!pr->betaA||!pr->betaB) die("OOM");

  int W=pm->win; if(W%2==0) W+=1; int h=(W-1)/2;
  for(int j=1;j<=pr->J;++j){
    int Lw=j-h; if(Lw<1) Lw=1;
    int Rw=j+h; if(Rw>pr->J) Rw=pr->J;
    if(Rw>L) Rw=L;
    if(Lw>Rw){ pr->betaA[j]=0; pr->betaB[j]=0; continue; }
    int sA=MA[Rw]-MA[Lw-1];
    int sB=MB[Rw]-MB[Lw-1];
    int s=sA-sB;
    pr->betaA[j] = - pm->kappa * s;
    pr->betaB[j] = + pm->kappa * s;
  }

  free(Au); free(Bu); free(MA); free(MB);
}

static void build_bands(Pre*pr, const Params*pm){
  pr->tHatA=(int*)malloc(sizeof(int)*(pr->J+1));
  pr->tHatB=(int*)malloc(sizeof(int)*(pr->J+1));
  pr->loA=(int*)malloc(sizeof(int)*(pr->J+1));
  pr->hiA=(int*)malloc(sizeof(int)*(pr->J+1));
  pr->loB=(int*)malloc(sizeof(int)*(pr->J+1));
  pr->hiB=(int*)malloc(sizeof(int)*(pr->J+1));
  if(!pr->tHatA||!pr->tHatB||!pr->loA||!pr->hiA||!pr->loB||!pr->hiB) die("OOM");

  int B = (pm->band>0)? pm->band : imax(10, pm->win);
  for(int j=1;j<=pr->J;++j){
    int ta = (j<=pr->La) ? pr->mapA[j] : iround( (double)j * pr->T / pr->J );
    int tb = (j<=pr->Lb) ? pr->mapB[j] : iround( (double)j * pr->T / pr->J );
    pr->tHatA[j]=ta; pr->tHatB[j]=tb;
  }
  int W=pm->win; if(W%2==0) W+=1; int h=(W-1)/2;
  int *tmpA=(int*)malloc(sizeof(int)*(pr->J+1));
  int *tmpB=(int*)malloc(sizeof(int)*(pr->J+1));
  if(!tmpA||!tmpB) die("OOM");
  for(int j=1;j<=pr->J;++j){
    int Lw=j-h; if(Lw<1) Lw=1;
    int Rw=j+h; if(Rw>pr->J) Rw=pr->J;
    long sumA=0,sumB=0; int cnt=0;
    for(int q=Lw;q<=Rw;++q){ sumA+=pr->tHatA[q]; sumB+=pr->tHatB[q]; cnt++; }
    tmpA[j]=(int)(sumA/cnt);
    tmpB[j]=(int)(sumB/cnt);
  }
  for(int j=1;j<=pr->J;++j){ pr->tHatA[j]=tmpA[j]; pr->tHatB[j]=tmpB[j]; }
  free(tmpA); free(tmpB);

  pr->loA[0]=0; pr->hiA[0]=pr->T; pr->loB[0]=0; pr->hiB[0]=pr->T;
  pr->loA[pr->J]=0; pr->hiA[pr->J]=pr->T; pr->loB[pr->J]=0; pr->hiB[pr->J]=pr->T;
  for(int j=1;j<pr->J;++j){
    int lo,hi;
    lo = pr->tHatA[j]-B; if(lo<0) lo=0; hi=pr->tHatA[j]+B; if(hi>pr->T) hi=pr->T;
    pr->loA[j]=lo; pr->hiA[j]=hi;
    lo = pr->tHatB[j]-B; if(lo<0) lo=0; hi=pr->tHatB[j]+B; if(hi>pr->T) hi=pr->T;
    pr->loB[j]=lo; pr->hiB[j]=hi;
  }
}

static inline int beta_at(const Pre*pr,int m,int j){ return (m==0)? pr->betaA[j] : pr->betaB[j]; }
static inline char parent_char(const Pre*pr,int t,int m){ return (m==0)? pr->A[t] : pr->B[t]; }
static inline double switch_cost_at(const Pre*pr,const Params*pm,int t){
  return (double)pm->sw_base + (pr->nearGap[t]? pm->gap_swell : 0);
}
static inline double sub_cost_at(const Pre*pr,const Params*pm,int t,int j,int m,int stay){
  char p = parent_char(pr,t,m);
  int mu = pr->mu[t];
  int match = (j>=1 && j<=pr->J && is_base(p) && pr->Y[j]==p);
  double cost = match? 0.0 : (double)mu;
  if(pr->I[t] && match && pm->reward>0) cost -= (double)pm->reward;
  if(stay && pm->stick>0) cost -= (double)pm->stick;
  cost += (double)beta_at(pr,m,j);
  return cost;
}
static inline double ins_cost_at(const Pre*pr,const Params*pm,int j,int m){
  return (double)pm->gap + (double)beta_at(pr,m,j);
}
static inline double del_cost_at(const Params*pm){ return (double)pm->gap; }

static inline int node_id(int T,int J,int t,int j,int m){ return ((t*(J+1)+j)<<1) | (m&1); }
static inline void node_ijm(int id,int J,int*pt,int*pj,int*pm){
  *pm = id & 1; int v = id>>1; *pj = v % (J+1); *pt = v / (J+1);
}

typedef struct { double *db; int *best_to; } BackDP;

static inline int allowed_t_for(int j,int m,const Pre*pr,int *plo,int *phi){
  if(j==0 || j==pr->J){ *plo=0; *phi=pr->T; return 1; }
  if(m==0){ *plo=pr->loA[j]; *phi=pr->hiA[j]; return (*plo<=*phi); }
  else    { *plo=pr->loB[j]; *phi=pr->hiB[j]; return (*plo<=*phi); }
}

static BackDP run_backward_dp_banded(const Pre*pr, const Params*pm){
  int T=pr->T, J=pr->J;
  int N = (T+1)*(J+1)*2;
  BackDP B;
  B.db = (double*)malloc(sizeof(double)*N);
  B.best_to = (int*)malloc(sizeof(int)*N);
  if(!B.db||!B.best_to) die("OOM");
  const double INF = 1e50;
  for(int i=0;i<N;++i){ B.db[i]=INF; B.best_to[i]=-1; }

  B.db[node_id(T,J,T,J,0)] = 0.0;
  B.db[node_id(T,J,T,J,1)] = 0.0;

  for(int j=J; j>=0; --j){
    for(int m=0;m<2;++m){
      int lo,hi;
      if(!allowed_t_for(j,m,pr,&lo,&hi)) continue;
      for(int t=hi; t>=lo; --t){
        if(t==T && j==J) continue;
        int id = node_id(T,J,t,j,m);

        double best = B.db[id];
        int bestto=-1;

        if(t<T && j<J){
          char p = parent_char(pr,t+1,m);
          if(is_base(p)){
            int to = node_id(T,J,t+1,j+1,m);
            double cand = sub_cost_at(pr,pm,t+1,j+1,m,1) + B.db[to];
            if(cand < best){ best=cand; bestto=to; }
          }
        }
        if(t<T && j<J){
          int mb=1-m;
          char p = parent_char(pr,t+1,mb);
          if(is_base(p)){
            int to = node_id(T,J,t+1,j+1,mb);
            double cand = switch_cost_at(pr,pm,t+1) + sub_cost_at(pr,pm,t+1,j+1,mb,0) + B.db[to];
            if(cand < best){ best=cand; bestto=to; }
          }
        }
        if(j<J){
          int to = node_id(T,J,t,j+1,m);
          double cand = ins_cost_at(pr,pm,j+1,m) + B.db[to];
          if(cand < best){ best=cand; bestto=to; }
        }
        if(t<T){
          int to = node_id(T,J,t+1,j,m);
          double cand = del_cost_at(pm) + B.db[to];
          if(cand < best){ best=cand; bestto=to; }
        }

        B.db[id]=best; B.best_to[id]=bestto;
      }
    }
  }
  return B;
}

typedef struct {
  uint8_t *orig; // 1..J (0=A,1=B)
  int switches, stays, insertions, deletions;
  int matches, mismatches;
  double total_cost;
} PathOut;

static void traceback_and_score(const BackDP*B, const Pre*pr, const Params*pm, PathOut*out){
  int T=pr->T, J=pr->J;
  int idA=node_id(T,J,0,0,0), idB=node_id(T,J,0,0,1);
  int start = (B->db[idA] <= B->db[idB]) ? idA : idB;

  memset(out->orig, 0, (size_t)J+1);
  out->switches=0; out->stays=0; out->insertions=0; out->deletions=0;
  out->matches=0; out->mismatches=0;
  out->total_cost=0.0;

  int cur=start;
  int t,j,m; node_ijm(cur,J,&t,&j,&m);

  while(!(t==T && j==J)){
    int to = B->best_to[cur];
    if(to<0) break;
    int t2,j2,m2; node_ijm(to,J,&t2,&j2,&m2);

    if(t2==t+1 && j2==j+1){
      char pc = parent_char(pr,t+1,m2);
      int match = (is_base(pc) && j+1>=1 && j+1<=J && pr->Y[j+1]==pc);
      if(m2==m){
        out->stays++;
        out->total_cost += sub_cost_at(pr,pm,t2,j2,m2,1);
      }else{
        out->switches++;
        out->total_cost += switch_cost_at(pr,pm,t2) + sub_cost_at(pr,pm,t2,j2,m2,0);
      }
      if(match) out->matches++; else out->mismatches++;
      if(j2==j+1) out->orig[j2]=(uint8_t)m2;
    }else if(t2==t && j2==j+1){
      out->insertions++;
      out->total_cost += ins_cost_at(pr,pm,j2,m2);
      out->orig[j2]=(uint8_t)m2;
    }else if(t2==t+1 && j2==j){
      out->deletions++;
      out->total_cost += del_cost_at(pm);
    }
    cur=to; t=t2; j=j2; m=m2;
  }
}

static void write_outputs(const char*outprefix, const char*child_name,
                          const Pre*pr, const PathOut*po)
{
  char path_orig[1024], path_sum[1024];
  child_name = "";
  snprintf(path_orig,sizeof(path_orig), "%s.%s.origins.tsv", outprefix, child_name);
  snprintf(path_sum, sizeof(path_sum ), "%s.%s.summary.tsv", outprefix, child_name);

  FILE*f1=fopen(path_orig,"w");
  if(!f1) die("open %s: %s", path_orig, strerror(errno));
  fprintf(f1, "j\tparent\n");
  for(int j=1;j<=pr->J;++j){
    char c = (po->orig[j]==0)? 'A':'B';
    fprintf(f1, "%d\t%c\n", j, c);
  }
  fclose(f1);

  FILE*f2=fopen(path_sum,"w");
  if(!f2) die("open %s: %s", path_sum, strerror(errno));
  // required format: header + one row
  fprintf(f2, "total\tmismatches\tinsertions\tdeletions\tswitches\n");
  long total_int = llround(po->total_cost);
  fprintf(f2, "%ld\t%d\t%d\t%d\t%d\n",
          total_int, po->mismatches, po->insertions, po->deletions, po->switches);
  fclose(f2);
}

static int run_once_with_band(const char*parents_path, const char*child_path,
                              const char*outprefix, Params *pm){
  Parents P = read_parents_aln(parents_path);
  Fasta   Y = read_fasta_one(child_path);

  Pre pr; memset(&pr,0,sizeof(pr));
  pr.T=P.A.len; pr.J=Y.len; pr.A=P.A.seq; pr.B=P.B.seq; pr.Y=Y.seq;

  build_informative_and_neargap(&pr, pm);
  build_ungap_maps(&pr);
  build_unary_bias(&pr, pm);
  build_bands(&pr, pm);

  BackDP B = run_backward_dp_banded(&pr, pm);

  int idA=node_id(pr.T,pr.J,0,0,0), idB=node_id(pr.T,pr.J,0,0,1);
  int reachable = (B.db[idA]<1e49 || B.db[idB]<1e49);

  if(reachable){
    PathOut po; po.orig=(uint8_t*)malloc((size_t)pr.J+1); if(!po.orig) die("OOM");
    traceback_and_score(&B, &pr, pm, &po);
    write_outputs(outprefix, Y.name, &pr, &po);
    free(po.orig);
  }

  free(B.db); free(B.best_to);
  free(pr.I); free(pr.mu); free(pr.nearGap);
  free(pr.betaA); free(pr.betaB);
  free(pr.mapA); free(pr.mapB);
  free(pr.tHatA); free(pr.tHatB);
  free(pr.loA); free(pr.hiA); free(pr.loB); free(pr.hiB);
  free(P.A.seq); free(P.B.seq); free(Y.seq);

  return reachable? 0 : 1;
}

static void parse_args(int argc,char**argv,
                       char **parents_path, char **child_path, char **outprefix,
                       Params *pm)
{
  set_default_params(pm);
  *parents_path=NULL; *child_path=NULL; *outprefix=NULL;
  for(int i=1;i<argc;++i){
    if(strcmp(argv[i],"--parents-aln")==0 && i+1<argc) *parents_path=argv[++i];
    else if(strcmp(argv[i],"--child")==0 && i+1<argc) *child_path=argv[++i];
    else if(strcmp(argv[i],"--out")==0 && i+1<argc) *outprefix=argv[++i];

    else if(strcmp(argv[i],"--mismatch")==0 && i+1<argc){ // red.c compatibility
      int v = atoi(argv[++i]); pm->mismatch_eq=v; pm->mismatch_inf=v;
    }
    else if(strcmp(argv[i],"--mu-eq")==0 && i+1<argc)  pm->mismatch_eq=atoi(argv[++i]);
    else if(strcmp(argv[i],"--mu-inf")==0 && i+1<argc) pm->mismatch_inf=atoi(argv[++i]);

    else if(strcmp(argv[i],"--reward")==0 && i+1<argc) pm->reward=atoi(argv[++i]);
    else if(strcmp(argv[i],"--gap")==0 && i+1<argc)    pm->gap=atoi(argv[++i]);
    else if(strcmp(argv[i],"--switch")==0 && i+1<argc) pm->sw_base=atoi(argv[++i]);
    else if(strcmp(argv[i],"--gap-swell")==0 && i+1<argc) pm->gap_swell=atoi(argv[++i]);
    else if(strcmp(argv[i],"--near-gap")==0 && i+1<argc)  pm->gap_near=atoi(argv[++i]);
    else if(strcmp(argv[i],"--gap-near")==0 && i+1<argc)  pm->gap_near=atoi(argv[++i]); // alias

    else if(strcmp(argv[i],"--stick")==0 && i+1<argc)  pm->stick=atoi(argv[++i]);
    else if(strcmp(argv[i],"--win")==0 && i+1<argc)    pm->win=atoi(argv[++i]);
    else if(strcmp(argv[i],"--kappa")==0 && i+1<argc)  pm->kappa=atoi(argv[++i]);
    else if(strcmp(argv[i],"--band")==0 && i+1<argc)   pm->band=atoi(argv[++i]);
    else die("Unknown/invalid arg: %s", argv[i]);
  }
  if(!*parents_path || !*child_path || !*outprefix) die("Required: --parents-aln, --child, --out");
  if(pm->win%2==0) pm->win += 1;
  if(pm->band==0) pm->band = -1;
  if(pm->band<0)  pm->band = imax(10, pm->win);
}

int main(int argc,char**argv){
  char *parents_path,*child_path,*outprefix; Params pm;
  parse_args(argc,argv,&parents_path,&child_path,&outprefix,&pm);

  int rc = run_once_with_band(parents_path, child_path, outprefix, &pm);
  if(rc==0) return 0;

  fprintf(stderr,"[WARN] unreachable with band=%d; widening and retry...\n", pm.band);
  int saved_band = pm.band;
  for(int step=0; step<2; ++step){
    pm.band = saved_band * 2;
    int rc2 = run_once_with_band(parents_path, child_path, outprefix, &pm);
    if(rc2==0) return 0;
    fprintf(stderr,"[WARN] still unreachable with band=%d; widening...\n", pm.band);
  }
  fprintf(stderr,"[ERROR] DP path unreachable even after widening band.\n");
  return 2;
}

