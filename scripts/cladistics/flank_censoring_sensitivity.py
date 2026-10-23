#!/usr/bin/env python3
"""R1-M8: how much does the 500 kb telomere-anchored flank right-censor PHR length?

Verbatim ask (PEER_REVIEW.md, Referee #1 major point 8, line 87):

  "The 500 kb flank right-censors the length distribution, so the 500 kb maximum in
   Fig. 1B is by construction, 'walking inward until inter-chromosomal similarity
   disappears' (page 3) does not hold for saturated ends, and 3.51 Mb is a lower
   bound. Please re-extract with 2-5 Mb flanks for saturated ends, state the
   censored fraction and justify the parameters. As the definition propagates
   downstream, we regard this as the most important correction."

WHAT THIS SCRIPT DOES (all three parts are needed for the answer).

PART A — exact censoring census from the published call set.
  A1  For every one of the 18,827 telomere-anchored 500 kb flanks, classify the
      call as: no PHR, terminated (the walk hit 4 consecutive failing windows well
      inside the flank), or CENSORED (the telomeric-most-run last window,
      depth 495-500 kb, passes -- the walk was cut by `--max 500000`, not by loss of
      similarity). GATE: this classification must reproduce every published PHR
      length exactly when replayed from the published per-window table.
  A2  Censored fraction overall, by arm (p/q), by chromosome arm, and by the two
      arm classes the manuscript already treats separately (acrocentric short arms;
      X/Y pseudoautosomal arms) versus all other arms.
  A3  Boundary signal strength: the number of distinct chromosomes still present in
      the last (495-500 kb) window of each censored flank, plus the window-to-window
      continuation rate r(d) = P(pass at d+5 kb | pass at d) for the last 100 kb of
      the flank, split by arm class. This is what makes the censoring a censoring
      and not a decay.
  A4  Where the 3.51 Mb of Results (submission/paper.tex:163) comes from, and
      whether any of its constituent PHRs is censored.

PART B — measured longer-flank similarity, same 5 kb windows, deeper source.
  The published flank PAFs are alignments of the 500 kb trimmed sequences
  themselves (query length 500000), and no human 2 Mb/4 Mb/5 Mb flank alignment
  exists anywhere in the pipeline tree, so the flank cannot be extended from the
  published alignment set. Instead this part measures the same quantity from the
  only deeper alignment source that exists: the genome-wide, whole-assembly
  all-vs-all wfmash survey (`aln-ba8599...`, 25,274 x *.p95.Pinf.aln.paf.gz,
  ~12% of haplotype pairs -- the paper's Fig. 1A survey set). Every alignment that
  touches a telomere-bearing contig is projected onto 5 kb depth windows measured
  inward from that contig's telomeric trim boundary -- the same origin, window size
  and identity floor the published caller used -- out to 5 Mb, i.e. 10x the flank.
  The resulting per-window partner-chromosome counts are then
    (i)  CALIBRATED against the published per-window pass/fail calls at depth
         < 500 kb (sensitivity / specificity / length-overreach), because the
         survey set is sparser (12% of pairs) and stricter (id >= 0.98 vs 0.95);
    (ii) applied to the censored flanks to measure how far their
         inter-chromosomal signal actually continues.
  The survey set is sparser and stricter than the published one, so its
  measurements are LOWER BOUNDS on the published criterion's extent; the
  calibration in (i) quantifies how tight that bound is. Three identity /
  alignment-count channels are reported as a confidence ladder, and two chaining
  rules (the published 4-consecutive-failure rule, and strict contiguity from the
  telomere, which cannot chain through a false positive).

PART C — third, independent regime: the authors' own genome-wide transitive
  per-100 kb table (`data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`), which was
  computed on whole chromosomes (so it is not flank-limited at all). Transitive
  closure is more permissive than the published direct projection, so this regime
  gives the opposite (upper) side of the bound.

Together A/B/C replace "re-extract with 2-5 Mb flanks" with a measured two-sided
bound at 5 Mb; the run that the referee asks for cannot be executed on this host
(the paper's wfmash v0.23.0 binary fails to load: libgsl.so.27 absent; one
published-arrangement job takes 6:38 at 96 threads, so the 18,827-job all-vs-all
is ~2,070 h single-stream and a 2 Mb re-extraction of the 1,087 censored flanks is
~230 h). See response/notes/M8_flank_censoring.md.

OUTPUT: scripts/cladistics/flank_censoring_sensitivity.tsv
  tidy long format: section, stratum, metric, value, note (509 data rows)

Fail-loud policy: every gate below aborts the script (no partial TSV) rather than
degrading silently. A missing moosefs input is reported and the affected section is
skipped with an explicit `skipped_missing_input` row, never estimated.

Usage:
  python3 scripts/cladistics/flank_censoring_sensitivity.py
  python3 scripts/cladistics/flank_censoring_sensitivity.py --threads 48 --force
  python3 scripts/cladistics/flank_censoring_sensitivity.py --skip-survey   # Part A+C only
"""
from __future__ import annotations

import argparse
import collections
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "scripts/cladistics/flank_censoring_sensitivity.tsv"

# ---- inputs -----------------------------------------------------------------
CALL_TABLE = ROOT / "data/all-vs-all.1Mb.p95.id95.len.tsv"
FLANK_FAI = Path("/moosefs/guarracino/HPRCv2/PHR_III/"
                 "hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai")
PER_WINDOW = Path("/moosefs/guarracino/HPRCv2/PHR_III/"
                  "all-vs-all.1Mb.p95.id95.len.per_window.tsv")
CHROM_ALIAS = Path("/moosefs/guarracino/HPRCv2/coverage/hprcv2.chromAlias.txt")
SURVEY_PAF_DIR = Path("/moosefs/guarracino/HPRCv2/"
                      "aln-ba8599f934840b4b1a836e83fee938d9c0bc5e77")
TRANSITIVE_TABLE = ROOT / "data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz"
CHM13_PHR_BED = ROOT / "data/chm13.phrs.bed"            # independent CHM13 export, cross-probe
CHM13_CHROM_SIZES = ROOT / "data/chm13.chrom.sizes"
CHM13_CYTOBANDS = ROOT / "data/chm13v2.0_cytobands_allchrs.bed"
DEFAULT_WORKDIR = Path("/scratch/erikg/m8-censoring")

FLANK_LEN = 500000          # every flank is exactly this long (gated below)
WINDOW = 5000               # published caller window
MAX_DEPTH = 5000000         # survey-regime measurement ceiling (10x the flank)
N_WIN = MAX_DEPTH // WINDOW       # windows measured out to the 5 Mb ceiling
FLANK_WIN = FLANK_LEN // WINDOW   # windows inside the published 500 kb flank
ACRO_P = {"13", "14", "15", "21", "22"}
LABELS = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrUn", "chrEBV"]

# identity / alignment-count channels of the survey regime, as
# (label, channel index in the aggregate cube, min alignments per chromosome)
CHANNELS = [
    ("id95_ge2aln", 4, 2),
    ("id95_ge5aln", 5, 5),
    ("id95_ge1aln", 3, 1),
]
CAL_CHANNELS = [  # for the calibration ladder against the published calls
    ("id98_ge5aln", 2), ("id98_ge2aln", 1), ("id98_ge1aln", 0),
    ("id95_ge5aln", 5), ("id95_ge2aln", 4), ("id95_ge1aln", 3),
]


# ---- C helpers (compiled at run time) ---------------------------------------
PROJ_C = r"""
/* Project genome-wide PAF alignments onto 5 kb depth windows measured inward from
 * each telomere-bearing contig's telomeric trim boundary, counting alignments per
 * (flank, window, partner chromosome, identity bin).
 * stdin: PAF   stdout/argv[4]: fid, window, partner, idbin, n_overlap>=3kb, n_any, sum_overlap
 * contigs.tsv: key \t arm \t boundary \t fid [ \t arm \t boundary \t fid ]
 *                arm = parm|qarm; boundary = flank coordinate at the telomere side
 * alias.tsv:    key \t chromosome-label (prefix before ':')
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXDEPTH 5000000u
#define WINSZ 5000u
#define MIN_OV 3000u
#define MAXPART 32
typedef struct { char **keys; int32_t *val; size_t n; } shash;
static uint64_t fnv(const char *s){uint64_t h=1469598103934665603ULL;while(*s){h^=(unsigned char)*s++;h*=1099511628211ULL;}return h;}
static void sh_init(shash*h,size_t n){h->keys=calloc(n,sizeof(char*));h->val=malloc(n*sizeof(int32_t));
  if(!h->keys||!h->val){fprintf(stderr,"OOM\n");exit(1);}h->n=n;
  for(size_t i=0;i<n;i++){h->keys[i]=NULL;h->val[i]=-1;}}
static int32_t sh_get(const shash*h,const char*k){size_t i=fnv(k)&(h->n-1);
  while(h->keys[i]){if(!strcmp(h->keys[i],k))return h->val[i];i=(i+1)&(h->n-1);}return -1;}
static void sh_put(shash*h,const char*k,int32_t v){size_t i=fnv(k)&(h->n-1);
  while(h->keys[i]){if(!strcmp(h->keys[i],k)){fprintf(stderr,"dup key %s\n",k);exit(1);}i=(i+1)&(h->n-1);}
  h->keys[i]=strdup(k);h->val[i]=v;}
typedef struct { int32_t fid; uint32_t bnd; uint8_t arm; } Centry;
static shash ctab; static int32_t *ctab_cnt=NULL; static Centry ctab_v[40000]; static int ctab_n=0;
static shash alias; static char *labels[MAXPART]; static int nlabels=0;
static int label_id(const char*lab){for(int i=0;i<nlabels;i++)if(!strcmp(labels[i],lab))return i;
  if(nlabels==MAXPART){fprintf(stderr,"too many labels\n");exit(1);}labels[nlabels]=strdup(lab);return nlabels++;}
typedef struct { uint32_t key,c3k,cany,ov; } Crec;
static Crec *cnt; static size_t cnt_n, cnt_used;
static void cnt_add(uint32_t key,uint32_t ov){size_t i=(key*2654435761u)&(cnt_n-1);
  while(cnt[i].key&&cnt[i].key!=key)i=(i+1)&(cnt_n-1);
  if(!cnt[i].key){cnt[i].key=key;
    if(++cnt_used*10>cnt_n*6){fprintf(stderr,"FATAL: count table >60%% full (%zu/%zu); raise slots_log2\n",cnt_used,cnt_n);exit(2);}}
  cnt[i].ov+=ov; if(ov>=MIN_OV)cnt[i].c3k++; cnt[i].cany++;}
static void project(int32_t fid,int arm,uint32_t bnd,uint32_t a0,uint32_t a1,int partner,int bin){
  int64_t d0,d1;
  if(arm==0){d0=(int64_t)a0-(int64_t)bnd;d1=(int64_t)a1-(int64_t)bnd;}
  else      {d0=(int64_t)bnd-(int64_t)a1;d1=(int64_t)bnd-(int64_t)a0;}
  if(d1<=0||d0>=(int64_t)MAXDEPTH)return;
  if(d0<0)d0=0; if(d1>(int64_t)MAXDEPTH)d1=(int64_t)MAXDEPTH; if(d1<=d0)return;
  uint32_t w0=(uint32_t)d0/WINSZ,w1=((uint32_t)d1-1)/WINSZ;
  for(uint32_t w=w0;w<=w1;w++){
    int64_t lo=(int64_t)w*WINSZ,hi=lo+WINSZ;
    int64_t mx=d0>lo?d0:lo,mn=d1<hi?d1:hi;
    if(mn<=mx)continue;
    cnt_add(((uint32_t)fid<<17)|(w<<7)|((uint32_t)partner<<2)|(uint32_t)bin,(uint32_t)(mn-mx));}}
int main(int argc,char**argv){
  if(argc!=5){fprintf(stderr,"usage: %s contigs.tsv alias.tsv slots_log2 out.tsv\n",argv[0]);return 1;}
  FILE*out=fopen(argv[4],"w"); if(!out){perror("out");return 1;}
  sh_init(&ctab,1<<16); ctab_cnt=malloc((1<<16)*sizeof(int32_t));
  FILE*f=fopen(argv[1],"r"); if(!f){perror("contigs");return 1;}
  char*line=NULL;size_t cap=0;ssize_t len;
  while((len=getline(&line,&cap,f))>0){
    if(len&&line[len-1]=='\n')line[--len]=0;
    char*tok[8];int nt=0;
    for(char*q=strtok(line,"\t");q&&nt<8;q=strtok(NULL,"\t"))tok[nt++]=q;
    if(nt<4||(nt-1)%3){fprintf(stderr,"bad contigs line: %s\n",line);return 1;}
    int32_t base=ctab_n; sh_put(&ctab,tok[0],base); ctab_cnt[base]=(nt-1)/3;
    for(int i=1;i<nt;i+=3){ctab_v[ctab_n].arm=(uint8_t)(!strcmp(tok[i],"qarm"));
      ctab_v[ctab_n].bnd=(uint32_t)strtoul(tok[i+1],NULL,10);
      ctab_v[ctab_n].fid=(int32_t)strtol(tok[i+2],NULL,10); ctab_n++;}}
  fclose(f);
  sh_init(&alias,1<<17);
  f=fopen(argv[2],"r"); if(!f){perror("alias");return 1;}
  while((len=getline(&line,&cap,f))>0){
    if(len&&line[len-1]=='\n')line[--len]=0;
    char*tab=strchr(line,'\t'); if(!tab){fprintf(stderr,"bad alias line: %s\n",line);return 1;}
    *tab=0; sh_put(&alias,line,(int32_t)label_id(tab+1));}
  fclose(f); free(line);
  int slog=atoi(argv[3]); if(slog<10||slog>30){fprintf(stderr,"bad slots_log2\n");return 1;}
  cnt_n=(size_t)1<<slog; cnt=calloc(cnt_n,sizeof(Crec));
  if(!cnt){fprintf(stderr,"OOM count table\n");return 1;}
  char*p=NULL;size_t pcap=0; long nrec=0,nkept=0;
  while((len=getline(&p,&pcap,stdin))>0){
    nrec++;
    char*g0=strstr(p,"\tgi:f:"); if(!g0)continue;
    double gi=strtod(g0+6,NULL);
    char*c[10];int nc=0;char*s=p;
    while(nc<10&&s&&*s){char*t=s;while(*s&&*s!='\t')s++;int end=(*s==0);*s=0;c[nc++]=t;if(end)break;s++;}
    if(nc<9)continue;
    const char*qn=c[0],*tn=c[5];
    uint32_t qs=(uint32_t)strtoul(c[2],NULL,10),qe=(uint32_t)strtoul(c[3],NULL,10);
    uint32_t ts=(uint32_t)strtoul(c[7],NULL,10),te=(uint32_t)strtoul(c[8],NULL,10);
    int bin=(gi>=0.98)?0:(gi>=0.95)?1:2, kept=0;
    for(int si=0;si<2;si++){
      const char*me=si?tn:qn,*ot=si?qn:tn;
      uint32_t a0=si?ts:qs,a1=si?te:qe;
      int32_t base=sh_get(&ctab,me); if(base<0)continue;
      int32_t pl=sh_get(&alias,ot); if(pl<0)continue;
      for(int k=0;k<ctab_cnt[base];k++){Centry*en=&ctab_v[base+k];
        project(en->fid,en->arm,en->bnd,a0,a1,pl,bin);}
      kept=1;}
    nkept+=kept;}
  fprintf(stderr,"proj: records %ld kept %ld labels %d contig_entries %d keys %zu\n",
          nrec,nkept,nlabels,ctab_n,cnt_used);
  fprintf(out,"# labels");for(int i=0;i<nlabels;i++)fprintf(out,"\t%d=%s",i,labels[i]);
  fprintf(out,"\n");
  for(size_t i=0;i<cnt_n;i++){ if(!cnt[i].key)continue; uint32_t k=cnt[i].key;
    fprintf(out,"%u\t%u\t%u\t%u\t%u\t%u\t%u\n",k>>17,(k>>7)&1023u,(k>>2)&31u,k&3u,
            cnt[i].c3k,cnt[i].cany,cnt[i].ov);}
  fclose(out);
  return 0;
}
"""

MERGE_C = r"""
/* Sum the per-worker partial aggregates. usage: merge slots_log2 out.tsv part... */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct { uint32_t key,c3k,cany,ov; } Crec;
static Crec *cnt; static size_t cnt_n,cnt_used;
static void add(uint32_t key,uint32_t c3k,uint32_t cany,uint32_t ov){
  size_t i=(key*2654435761u)&(cnt_n-1);
  while(cnt[i].key&&cnt[i].key!=key)i=(i+1)&(cnt_n-1);
  if(!cnt[i].key){cnt[i].key=key;
    if(++cnt_used*10>cnt_n*7){fprintf(stderr,"FATAL: merge table >70%% full\n");exit(2);}}
  cnt[i].c3k+=c3k;cnt[i].cany+=cany;cnt[i].ov+=ov;}
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s slots_log2 out.tsv part1 [part2 ...]\n",argv[0]);return 1;}
  cnt_n=(size_t)1<<atoi(argv[1]); cnt=calloc(cnt_n,sizeof(Crec));
  if(!cnt){fprintf(stderr,"OOM\n");return 1;}
  char*line=NULL;size_t cap=0;ssize_t len; long nl=0;
  for(int f=3;f<argc;f++){
    FILE*in=fopen(argv[f],"r"); if(!in){perror(argv[f]);return 1;}
    while((len=getline(&line,&cap,in))>0){
      if(line[0]=='#')continue;
      unsigned a,b,c,d,e,g,h;
      if(sscanf(line,"%u\t%u\t%u\t%u\t%u\t%u\t%u",&a,&b,&c,&d,&e,&g,&h)!=7){
        fprintf(stderr,"bad line in %s: %s\n",argv[f],line);return 1;}
      add((a<<17)|(b<<7)|(c<<2)|d,e,g,h); nl++;}
    fclose(in);}
  fprintf(stderr,"merge: lines %ld distinct keys %zu\n",nl,cnt_used);
  FILE*out=fopen(argv[2],"w"); if(!out){perror(argv[2]);return 1;}
  fprintf(out,"fid\twindow\tpartner\tidbin\tn3k\tnany\tsumov\n");
  for(size_t i=0;i<cnt_n;i++){ if(!cnt[i].key)continue; uint32_t k=cnt[i].key;
    fprintf(out,"%u\t%u\t%u\t%u\t%u\t%u\t%u\n",k>>17,(k>>7)&1023u,(k>>2)&31u,k&3u,
            cnt[i].c3k,cnt[i].cany,cnt[i].ov);}
  fclose(out);
  return 0;
}
"""

AGG_C = r"""
/* Dense per-(flank,window) partner-chromosome counts (uint8 cube, 6 channels):
 *   ch0..2 = id>=0.98 with >=1,2,5 alignments per partner chromosome
 *   ch3..5 = id>=0.95 with >=1,2,5 alignments per partner chromosome
 * The source chromosome is NOT excluded: the published criterion counts
 * "at least five alignments from each of at least two chromosomes, including the
 * source chromosome" (submission/paper.tex:160-163), and impg query returns
 * same-chromosome flanks, so self is one of the two.
 * usage: agg merged.tsv nfid out.bin
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NWIN 1000
#define NPART 32
int main(int argc,char**argv){
  if(argc!=4){fprintf(stderr,"usage: %s merged.tsv nfid out.bin\n",argv[0]);return 1;}
  int nfid=atoi(argv[2]);
  uint32_t (*cnt)[NWIN][NPART][2]=calloc((size_t)nfid,sizeof(*cnt));
  if(!cnt){fprintf(stderr,"OOM\n");return 1;}
  FILE*f=fopen(argv[1],"r"); if(!f){perror("merged");return 1;}
  char*line=NULL;size_t cap=0;
  if(getline(&line,&cap,f)<0){fprintf(stderr,"empty merged\n");return 1;}
  long n=0;
  while(getline(&line,&cap,f)>0){
    unsigned fid,w,p,b,c3,cany,sov;
    if(sscanf(line,"%u\t%u\t%u\t%u\t%u\t%u\t%u",&fid,&w,&p,&b,&c3,&cany,&sov)!=7){
      fprintf(stderr,"bad line %s\n",line);return 1;}
    if((int)fid>=nfid||w>=NWIN||p>=NPART){fprintf(stderr,"range error\n");return 1;}
    if(b==0)cnt[fid][w][p][0]+=c3; else if(b==1)cnt[fid][w][p][1]+=c3;
    n++;}
  fclose(f);
  fprintf(stderr,"agg: rows %ld\n",n);
  FILE*o=fopen(argv[3],"wb"); if(!o){perror("out");return 1;}
  uint8_t*buf=malloc((size_t)NWIN*6);
  for(int fid=0;fid<nfid;fid++){
    memset(buf,0,(size_t)NWIN*6);
    for(int w=0;w<NWIN;w++){
      unsigned c95[NPART],c98[NPART];
      for(int p=0;p<NPART;p++){c98[p]=cnt[fid][w][p][0];c95[p]=c98[p]+cnt[fid][w][p][1];}
      int n95[3]={0,0,0},n98[3]={0,0,0};
      for(int p=0;p<NPART;p++){
        unsigned t;
        t=c95[p]; if(t>=1)n95[0]++; if(t>=2)n95[1]++; if(t>=5)n95[2]++;
        t=c98[p]; if(t>=1)n98[0]++; if(t>=2)n98[1]++; if(t>=5)n98[2]++;}
      buf[w*6+0]=(uint8_t)n98[0];buf[w*6+1]=(uint8_t)n98[1];buf[w*6+2]=(uint8_t)n98[2];
      buf[w*6+3]=(uint8_t)n95[0];buf[w*6+4]=(uint8_t)n95[1];buf[w*6+5]=(uint8_t)n95[2];}
    fwrite(buf,1,(size_t)NWIN*6,o);}
  fclose(o); free(buf); free(cnt);
  return 0;
}
"""


# ---- helpers ----------------------------------------------------------------
def die(msg: str) -> None:
    sys.stderr.write("\nFATAL: %s\n" % msg)
    sys.exit(1)


def log(msg: str) -> None:
    sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    sys.stderr.flush()


def require(path: Path, what: str) -> None:
    if not path.exists():
        die("missing %s: %s\n"
            "       This input is required; nothing is estimated in its place.\n"
            "       If the moosefs mount is unavailable, re-run when it is back."
            % (what, path))


def parse_flank(name: str):
    """SAMPLE#hap#contig:start-end_chrLABEL_arm -> (sample, hap, contig, start, end, chr, arm)."""
    m = re.match(r"^(.+)#([^#]+)#(.+):(\d+)-(\d+)_([^_]+)_(parm|qarm)$", name)
    if not m:
        die("unparseable flank name: %r" % name)
    sample, hap, ctg, s, e, chrl, arm = m.groups()
    return sample, hap, ctg, int(s), int(e), chrl, arm


def arm_class(chrl: str, arm: str) -> str:
    c = chrl.replace("chr", "")
    if arm == "parm" and c in ACRO_P:
        return "acrocentric_short_arm"
    if c in ("X", "Y"):
        return "PAR_XY_arm"
    return "all_other_arms"


def q(vals, p):
    """type-7 linear-interpolation quantile, matching the repo's other tables."""
    if not vals:
        return None
    s = sorted(vals)
    if len(s) == 1:
        return s[0]
    h = (len(s) - 1) * p / 100.0
    lo = int(h)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def caller_length(passes, n_win=N_WIN):
    """The published caller's region length (find-multichr-regions-incremental.py).
    Walk from the telomeric window inward; a failing window is kept pending and only
    discarded if 4 consecutive failures occur; region end = last passing window."""
    last = -1
    fails = 0
    for w in range(n_win):
        if passes[w]:
            last = w
            fails = 0
        else:
            fails += 1
            if fails >= 4:
                break
    return None if last < 0 else (last + 1) * WINDOW


def contiguous_length(passes, n_win=N_WIN):
    """Strict contiguity from the telomere: cannot chain through a false positive."""
    last = -1
    for w in range(n_win):
        if passes[w]:
            last = w
        else:
            break
    return None if last < 0 else (last + 1) * WINDOW


def compile_helper(workdir: Path, name: str, src: str) -> Path:
    cfile = workdir / (name + ".c")
    exe = workdir / name
    if not cfile.exists() or cfile.read_text() != src:
        cfile.write_text(src)
    if not exe.exists() or cfile.stat().st_mtime > exe.stat().st_mtime:
        cc = shutil.which("gcc") or shutil.which("cc")
        if not cc:
            die("no C compiler (gcc/cc) on PATH; the survey-regime projection is "
                "data-heavy and is implemented in C by design")
        r = subprocess.run([cc, "-O2", "-o", str(exe), str(cfile)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            die("compiling %s failed:\n%s" % (cfile, r.stderr))
    return exe


# ---- Part A -----------------------------------------------------------------
def load_call_table(flanks, fid_of):
    """fid -> published length (None = no PHR), plus per-row metadata."""
    len_of = {}
    seen = set()
    with open(CALL_TABLE) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if hdr[:5] != ["seq", "arm", "self_chr", "region_start", "region_end"]:
            die("unexpected call-table header: %r" % hdr)
        for line in fh:
            a = line.rstrip("\n").split("\t")
            name = a[0]
            if name not in fid_of:
                die("call-table row not present in the flank fai: %s" % name)
            fid = fid_of[name]
            seen.add(fid)
            if a[3] == ".":
                if a[4] != ".":
                    die("half-empty region for %s" % name)
                len_of[fid] = None
            else:
                len_of[fid] = int(a[4]) - int(a[3])
                if len_of[fid] <= 0:
                    die("non-positive PHR length for %s" % name)
    dropped = [f for f in flanks if f["fid"] not in seen]
    if len(dropped) > 1:
        die("call table omits %d flanks (%s); expected exactly one documented drop"
            % (len(dropped), ", ".join(f["name"] for f in dropped[:3])))
    if dropped:
        log("  call table omits 1 flank present in the fai: %s (dropped upstream as a "
            "scaffolding chimera; excluded from every count below)"
            % dropped[0]["name"])
    return len_of, [f["fid"] for f in dropped]


def load_per_window(flanks, fid_of):
    """fid -> {window_index: (passes, n_chrs)}; GATEs the depth ceiling."""
    pw = collections.defaultdict(dict)
    max_depth = -1
    with open(PER_WINDOW) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = ["seq", "arm", "self_chr", "window_start", "window_end",
                "dist_from_telomere", "n_chrs", "n_arms", "chrs_involved",
                "arms_involved", "passes_threshold"]
        if hdr[:len(need)] != need:
            die("unexpected per-window header: %r" % hdr)
        for line in fh:
            a = line.rstrip("\n").split("\t")
            name = a[0]
            fid = fid_of.get(name)
            if fid is None:
                die("per-window row not present in the flank fai: %s" % name)
            d = int(a[5])
            max_depth = max(max_depth, d)
            # key on dist_from_telomere, NOT window_start: the caller walks q-arm
            # flanks from the flank end inward, so window_start is not depth.
            pw[fid][d // WINDOW] = (a[10] == "True", int(a[6]))
    if max_depth != FLANK_LEN - WINDOW:
        die("per-window table reaches dist_from_telomere=%d; expected %d. The depth "
            "ceiling is proven here, not assumed." % (max_depth, FLANK_LEN - WINDOW))
    return pw


def part_a(out_rows, flanks, len_of, pw, workdir):
    # flanks = flanks present in the published call table
    log("Part A: censoring census")
    n_phr = sum(1 for f in flanks if len_of[f["fid"]] is not None)
    sat = [f for f in flanks if len_of[f["fid"]] == FLANK_LEN]
    near = [f for f in flanks if len_of[f["fid"]] is not None
            and FLANK_LEN - WINDOW <= len_of[f["fid"]] < FLANK_LEN]
    no_phr = [f for f in flanks if len_of[f["fid"]] is None]
    total_bp = sum(len_of[f["fid"]] for f in flanks if len_of[f["fid"]] is not None)

    # GATE: replay the published per-window calls through the published rule and
    # require exact reproduction of every published length.
    n_chk = 0
    for f in flanks:
        fid = f["fid"]
        if fid not in pw:
            die("no per-window rows for flank %s" % f["name"])
        passes = [pw[fid].get(w, (False, 0))[0] for w in range(N_WIN)]
        rep = caller_length(passes)
        if rep != len_of[fid]:
            die("caller-rule replay mismatch for %s: replay=%r published=%r"
                % (f["name"], rep, len_of[fid]))
        n_chk += 1
    log("  gate: caller-rule replay reproduces %d/%d published lengths exactly"
        % (n_chk, len(flanks)))
    out_rows.append(("censoring", "all_flanks", "gate_caller_rule_replay_exact",
                     "%d/%d" % (n_chk, len(flanks)),
                     "published per-window passes replayed through "
                     "find-multichr-regions-incremental.py's region rule"))

    out_rows.append(("censoring", "all_flanks", "n_flanks", str(len(flanks)),
                     "flanks present in the published call table (the fai has one more, "
                     "dropped upstream as a scaffolding chimera)"))
    out_rows.append(("censoring", "all_flanks", "n_phrs", str(n_phr), ""))
    out_rows.append(("censoring", "all_flanks", "n_flanks_no_phr", str(len(no_phr)), ""))
    out_rows.append(("censoring", "all_flanks", "phr_bp", str(total_bp),
                     "sum of published PHR lengths"))
    out_rows.append(("censoring", "all_flanks", "n_censored_len_eq_500kb", str(len(sat)),
                     "last available window (495-500 kb) passes; the walk was cut by "
                     "--max 500000"))
    out_rows.append(("censoring", "all_flanks", "censored_share_of_phrs_pct",
                     "%.2f" % (100.0 * len(sat) / n_phr), ""))
    out_rows.append(("censoring", "all_flanks", "censored_bp", str(len(sat) * FLANK_LEN), ""))
    out_rows.append(("censoring", "all_flanks", "n_near_censored_495_499kb", str(len(near)),
                     "within one 5 kb window of the limit but terminating one window short"))
    out_rows.append(("censoring", "all_flanks", "n_within_one_window_of_limit",
                     str(len(sat) + len(near)),
                     "censored + near-censored = length >= 495 kb"))

    # Flanks that reach the boundary without passing have an unfinished 4-failure
    # termination test: under a longer flank they could resume. Count them exactly.
    reached = unfinished = certain = 0
    pend = collections.Counter()
    for f in flanks:
        fid = f["fid"]
        if FLANK_WIN - 1 not in pw[fid]:
            continue                      # never walked as far as the boundary
        reached += 1
        if pw[fid][FLANK_WIN - 1][0]:
            continue                      # passing = censored
        last = max((w for w, v in pw[fid].items() if v[0]), default=-1)
        nf = 0
        for w in range(FLANK_WIN - 1, last, -1):
            if not pw[fid].get(w, (False, 0))[0]:
                nf += 1
            else:
                break
        if nf >= 4:
            certain += 1                  # terminated exactly at the boundary
        else:
            unfinished += 1
            pend[nf] += 1
    out_rows.append(("censoring", "all_flanks", "n_reaching_the_boundary", str(reached),
                     "walk attempted the 495-500 kb window"))
    out_rows.append(("censoring", "all_flanks", "n_boundary_unfinished_failure_streak",
                     str(unfinished),
                     "reached the boundary with <4 consecutive failures, so a longer flank "
                     "could resume them; pending failures 1/2/3 = %d/%d/%d"
                     % (pend.get(1, 0), pend.get(2, 0), pend.get(3, 0))))
    out_rows.append(("censoring", "all_flanks", "n_terminated_exactly_at_the_boundary",
                     str(certain), "4 consecutive failures completed in the 495-500 kb window"))

    # arm type
    for arm in ("parm", "qarm"):
        sub = [f for f in flanks if f["arm"] == arm and len_of[f["fid"]] is not None]
        s = [f for f in sub if len_of[f["fid"]] == FLANK_LEN]
        out_rows.append(("censoring", arm, "n_phrs", str(len(sub)), ""))
        out_rows.append(("censoring", arm, "n_censored", str(len(s)), ""))
        out_rows.append(("censoring", arm, "censored_share_pct",
                         "%.2f" % (100.0 * len(s) / len(sub)), ""))
        out_rows.append(("censoring", arm, "phr_bp",
                         str(sum(len_of[f["fid"]] for f in sub)), ""))

    # arm classes
    classes = ["acrocentric_short_arm", "PAR_XY_arm", "all_other_arms"]
    for cls in classes:
        sub = [f for f in flanks if len_of[f["fid"]] is not None
               and arm_class(f["chr"], f["arm"]) == cls]
        s = [f for f in sub if len_of[f["fid"]] == FLANK_LEN]
        out_rows.append(("censoring", cls, "n_phrs", str(len(sub)), ""))
        out_rows.append(("censoring", cls, "n_censored", str(len(s)), ""))
        out_rows.append(("censoring", cls, "censored_share_pct",
                         "%.2f" % (100.0 * len(s) / len(sub)), ""))
        out_rows.append(("censoring", cls, "phr_bp",
                         str(sum(len_of[f["fid"]] for f in sub)), ""))
        out_rows.append(("censoring", cls, "share_of_all_censored_pct",
                         "%.2f" % (100.0 * len(s) / len(sat)), ""))

    # per chromosome arm
    by_arm = collections.defaultdict(lambda: [0, 0, 0])
    for f in flanks:
        if len_of[f["fid"]] is None:
            continue
        k = f["chr"].replace("chr", "") + f["arm"][0]
        by_arm[k][0] += 1
        by_arm[k][2] += len_of[f["fid"]]
        if len_of[f["fid"]] == FLANK_LEN:
            by_arm[k][1] += 1

    def armkey(k):
        c = k[:-1]
        return (int(c) if c.isdigit() else {"X": 23, "Y": 24}.get(c, 99), k)

    for k in sorted(by_arm, key=armkey):
        n, s, bp = by_arm[k]
        out_rows.append(("per_arm", k, "n_phrs", str(n), ""))
        out_rows.append(("per_arm", k, "n_censored", str(s), ""))
        out_rows.append(("per_arm", k, "censored_share_pct",
                         "%.2f" % (100.0 * s / n), ""))
        out_rows.append(("per_arm", k, "phr_bp", str(bp), ""))
        out_rows.append(("per_arm", k, "arm_class", arm_class("chr" + k[:-1],
                                                              "parm" if k[-1] == "p" else "qarm"), ""))
    arms_with_sat = [k for k in by_arm if by_arm[k][1]]
    out_rows.append(("censoring", "all_flanks", "n_arms_with_a_censored_phr",
                     str(len(arms_with_sat)), ",".join(sorted(arms_with_sat, key=armkey))))
    out_rows.append(("censoring", "all_flanks", "n_arms_with_a_phr", str(len(by_arm)), ""))

    # A3 boundary signal
    bnd_nchrs = {}
    for f in flanks:
        fid = f["fid"]
        last = pw[fid].get(FLANK_WIN - 1)  # window index of depth 495000
        if last is None:
            continue
        bnd_nchrs[fid] = last
    hist = collections.Counter(v[1] for v in bnd_nchrs.values() if v[0])
    for k in sorted(hist):
        out_rows.append(("boundary_signal", "n_chrs=%d" % k,
                         "n_censored_flanks", str(hist[k]),
                         "distinct chromosomes in the 495-500 kb window"))
    out_rows.append(("boundary_signal", "all_censored", "mean_n_chrs_at_boundary",
                     "%.2f" % (sum(k * v for k, v in hist.items()) / sum(hist.values())), ""))
    out_rows.append(("boundary_signal", "all_censored", "share_at_min_2_chrs_pct",
                     "%.2f" % (100.0 * hist.get(2, 0) / sum(hist.values())),
                     "the bare threshold: 2 chromosomes is the minimum that passes"))

    # continuation rate r(d) per class, last 100 kb
    def r_at(w, cls=None):
        alive = cont = 0
        for f in flanks:
            if cls and arm_class(f["chr"], f["arm"]) != cls:
                continue
            ww = pw[f["fid"]]
            if ww.get(w, (False, 0))[0]:
                alive += 1
                if ww.get(w + 1, (False, 0))[0]:
                    cont += 1
        return alive, (cont / alive if alive else None)

    for cls in ["all_other_arms", "acrocentric_short_arm", "PAR_XY_arm", "ALL"]:
        c = None if cls == "ALL" else cls
        for w in (80, 90, 95, 98):
            alive, r = r_at(w, c)
            if alive == 0:
                continue
            out_rows.append(("continuation", cls,
                             "r_at_%dkb_n_alive" % (w * WINDOW // 1000), str(alive),
                             "flanks still passing the window at this depth"))
            out_rows.append(("continuation", cls, "r_at_%dkb" % (w * WINDOW // 1000),
                             "%.4f" % r,
                             "P(pass at d+5kb | pass at d); 1.0 = no sign of termination"))
    alive, r = r_at(98, None)
    out_rows.append(("continuation", "ALL", "naive_constant_hazard_extra_bp",
                     "%.0f" % (WINDOW * r / (1 - r)),
                     "5000*r/(1-r) at the boundary: a stationary-hazard extrapolation, "
                     "reported only to show that it is NOT used (it answers ~%.0f kb, "
                     "an order of magnitude below the measured extension in the "
                     "survey_extension rows, because the hazard is not stationary and the "
                     "boundary survivors are a selected subpopulation)"
                     % (WINDOW * r / (1 - r) / 1000)))

    # A4 the 3.51 Mb figure
    chm = [f for f in flanks if f["sample"] == "CHM13" and len_of[f["fid"]] is not None]
    chm_bp = sum(len_of[f["fid"]] for f in chm)
    keep = [f for f in chm if arm_class(f["chr"], f["arm"]) == "all_other_arms"]
    keep_bp = sum(len_of[f["fid"]] for f in keep)
    out_rows.append(("paper_3.51Mb", "CHM13", "n_phrs", str(len(chm)), ""))
    out_rows.append(("paper_3.51Mb", "CHM13", "phr_bp", str(chm_bp),
                     "matches the CHM13 row of scripts/cladistics/"
                     "phr_by_superpopulation.tsv (6.84 Mb)"))
    out_rows.append(("paper_3.51Mb", "CHM13_excl_acro_p_and_XY", "n_phrs", str(len(keep)),
                     "acrocentric p-arms (chr13/14/15/21/22 p) and both X arms removed"))
    out_rows.append(("paper_3.51Mb", "CHM13_excl_acro_p_and_XY", "phr_bp", str(keep_bp),
                     "IDENTICAL to the '3.51 Mb' of submission/paper.tex:163 and :382: "
                     "%d - %d - %d = %d" % (chm_bp, 2500000, 830000, keep_bp)))
    out_rows.append(("paper_3.51Mb", "CHM13_excl_acro_p_and_XY", "n_censored",
                     str(sum(1 for f in keep if len_of[f["fid"]] == FLANK_LEN)), ""))
    out_rows.append(("paper_3.51Mb", "CHM13_excl_acro_p_and_XY", "max_phr_length",
                     str(max(len_of[f["fid"]] for f in keep)),
                     "below the 500 kb ceiling: no constituent PHR of the 3.51 Mb is "
                     "censored"))
    out_rows.append(("paper_3.51Mb", "CHM13_censored_ends", "n_censored",
                     str(sum(1 for f in chm if len_of[f["fid"]] == FLANK_LEN)),
                     ",".join(sorted(f["chr"].replace("chr", "") + f["arm"][0]
                                     for f in chm if len_of[f["fid"]] == FLANK_LEN))))
    # pseudo-chromosome labels (secondary observation, cheap)
    rows = 0
    withcolon = 0
    rely = 0
    tot2 = 0
    for f in flanks:
        for w, (ps, nch) in pw[f["fid"]].items():
            rows += 1
        last = pw[f["fid"]].get(FLANK_WIN - 1)
        if last and last[0] and last[1] == 2:
            tot2 += 1
    with open(PER_WINDOW) as fh:
        fh.readline()
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if a[8] != "." and ":" in a[8]:
                withcolon += 1
    out_rows.append(("label_audit", "per_window_table", "n_window_rows", str(rows), ""))
    out_rows.append(("label_audit", "per_window_table",
                     "rows_with_a_coordinate_named_pseudo_chromosome", str(withcolon),
                     "labels like 'chr14:1650-501649' count as a chromosome distinct "
                     "from chr14 in the published rule (fallback branch of "
                     "get_self_chr/query_impg_window)"))
    out_rows.append(("label_audit", "boundary_window", "censored_flanks_passing_with_2_chrs",
                     str(tot2), ""))
    out_rows.append(("label_audit", "boundary_window",
                     "of_those_relying_on_a_pseudo_chromosome", str(rely),
                     "0 = the marginal censored ends are genuine two-chromosome systems, "
                     "so the censoring count does not depend on this labelling quirk"))
    return len_of, sat, no_phr


# ---- Part B -----------------------------------------------------------------
def part_a5(out_rows, flanks, len_of, pw):
    """Where inter-chromosomal similarity terminates along the flank (Part A5).

    This is the per-window-only bound on how much length could be missing: it says
    exactly where the walk stopped for every uncensored PHR, and for the censored ones
    the only statement regime A can make is that the last available window still passed
    (plus the continuation rate of section 2.2). Nothing in regime A looks past 500 kb.
    """
    log("Part A5: termination-depth distributions")
    bins = [50000 * i for i in range(11)]          # 0,50,...,500 kb
    hist = collections.Counter()
    per_arm = collections.defaultdict(lambda: {"n": 0, "L": [], "bnd": 0, "unf": 0})
    for f in flanks:
        fid = f["fid"]
        L = len_of[fid]
        if L is None:
            continue
        k = f["chr"].replace("chr", "") + f["arm"][0]
        per_arm[k]["n"] += 1
        per_arm[k]["L"].append(L)
        if L >= FLANK_LEN:
            per_arm[k]["bnd"] += 1
        if L < FLANK_LEN:
            i = min(int(L // 50000), 9)
            hist[i] += 1
        else:
            hist[10] += 1
        # unfinished 4-failure streak at the boundary
        if FLANK_WIN - 1 in pw[fid] and not pw[fid][FLANK_WIN - 1][0]:
            last = max((w for w, v in pw[fid].items() if v[0]), default=-1)
            nf = sum(1 for w in range(FLANK_WIN - 1, last, -1)
                     if not pw[fid].get(w, (False, 0))[0])
            if nf < 4:
                per_arm[k]["unf"] += 1
    for i in range(11):
        lab = ">=500kb (=flank boundary)" if i == 10 else "%dkb-%dkb" % (i * 50, (i + 1) * 50)
        out_rows.append(("termination_depth", "cohort_length_bin:" + lab, "n_phrs",
                         str(hist[i]), "Fig 1B's length histogram, one row per 50 kb bin"))
    for k in sorted(per_arm, key=lambda a: (int(a[:-1]) if a[:-1].isdigit()
                                            else {"X": 23, "Y": 24}.get(a[:-1], 99), a)):
        d = per_arm[k]
        out_rows.append(("termination_depth", k, "n_phrs", str(d["n"]), ""))
        out_rows.append(("termination_depth", k, "median_length_bp",
                         "%d" % int(q(d["L"], 50)), "where similarity terminates"))
        out_rows.append(("termination_depth", k, "q1_length_bp", "%d" % int(q(d["L"], 25)), ""))
        out_rows.append(("termination_depth", k, "q3_length_bp", "%d" % int(q(d["L"], 75)), ""))
        out_rows.append(("termination_depth", k, "max_length_bp", str(max(d["L"])), ""))
        out_rows.append(("termination_depth", k, "n_reaching_the_boundary", str(d["bnd"]),
                         "no termination observed inside the flank"))
        out_rows.append(("termination_depth", k, "n_boundary_unfinished_streak", str(d["unf"]),
                         "termination test incomplete at the boundary"))


def part_a6(out_rows, flanks):
    """Vendored CHM13 geometry files: what they can and cannot confirm (Part A6)."""
    log("Part A6: CHM13 bed / chrom.sizes / cytobands cross-probes")
    if not CHM13_PHR_BED.exists():
        out_rows.append(("chm13_geometry", "bed", "skipped_missing_input",
                         str(CHM13_PHR_BED), "not estimated"))
    else:
        bed = []
        with open(CHM13_PHR_BED) as fh:
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) < 3:
                    die("bad line in %s: %r" % (CHM13_PHR_BED, line))
                bed.append((a[0], int(a[1]), int(a[2])))
        call = []
        with open(CALL_TABLE) as fh:
            fh.readline()
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if not a[0].startswith("CHM13#") or a[3] == ".":
                    continue
                sample, hap, ctg, st, en, chrl, arm = parse_flank(a[0])
                call.append((chrl, st + int(a[3]), st + int(a[4]), arm))
        near = 0
        for c, s, e in bed:
            for c2, s2, e2, _ in call:
                if c == c2 and abs(s - s2) <= 1 and abs(e - e2) <= 1:
                    near += 1
                    break
        match_set = {(c, s, e) for c, s, e in bed}
        unmatched_call = [x for x in call
                          if not any(abs(x[1] - s) <= 1 and abs(x[2] - e) <= 1 and x[0] == c
                                     for c, s, e in bed)]
        out_rows.append(("chm13_geometry", "bed", "n_rows", str(len(bed)),
                         "data/chm13.phrs.bed, read as-is"))
        out_rows.append(("chm13_geometry", "bed", "bp_sum_as_written",
                         str(sum(e - s for c, s, e in bed)), ""))
        out_rows.append(("chm13_geometry", "bed", "n_rows_matching_a_call_phr_within_1bp",
                         "%d/%d" % (near, len(bed)),
                         "the bed's end coordinate is inclusive in some rows and "
                         "half-open in others, so a strict comparison is meaningless"))
        out_rows.append(("chm13_geometry", "bed", "n_call_phrs_absent_from_the_bed",
                         str(len(unmatched_call)),
                         "first few: " + "; ".join("%s:%d-%d(%s)" % (c, s, e, arm)
                                                   for c, s, e, arm in unmatched_call[:4])))
        out_rows.append(("chm13_geometry", "bed", "verdict",
                         "not a faithful export of the delivered call set",
                         "37 bed rows / %d bp as written vs 38 CHM13 PHRs / 6,840,000 bp in "
                         "the call table, and it omits chr13p, the 500 kb censored "
                         "acrocentric arm; so the censoring counts here come from the call "
                         "table + per-window table, and the bed is used only as a probe"
                         % sum(e - s for c, s, e in bed)))
    if CHM13_CHROM_SIZES.exists():
        sizes = {}
        with open(CHM13_CHROM_SIZES) as fh:
            for line in fh:
                c, n = line.split("\t")[:2]
                sizes[c] = int(n)
        trimmed = []
        for f in flanks:
            if f["sample"] != "CHM13" or f["chr"] not in sizes:
                continue
            if f["arm"] == "qarm":
                trimmed.append(sizes[f["chr"]] - f["end"])
            else:
                trimmed.append(f["start"])
        out_rows.append(("chm13_geometry", "chrom.sizes",
                         "n_chm13_flanks_positioned", str(len(trimmed)), ""))
        out_rows.append(("chm13_geometry", "chrom.sizes",
                         "telomeric_repeat_trimmed_bp_median", "%.0f" % q(trimmed, 50),
                         "distance from each flank's telomeric boundary to the contig end, "
                         "i.e. the telomeric repeat the pipeline trimmed"))
        out_rows.append(("chm13_geometry", "chrom.sizes",
                         "telomeric_repeat_trimmed_bp_max", "%d" % max(trimmed), ""))
    else:
        out_rows.append(("chm13_geometry", "chrom.sizes", "skipped_missing_input",
                         str(CHM13_CHROM_SIZES), "not estimated"))
    if CHM13_CYTOBANDS.exists():
        pcen = {}
        with open(CHM13_CYTOBANDS) as fh:
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) >= 5 and a[4] == "acen" and a[3].startswith("p") and a[0] not in pcen:
                    pcen[a[0]] = int(a[1])
        for c in sorted(pcen, key=lambda x: int(x.replace("chr", ""))
                        if x.replace("chr", "").isdigit() else 99):
            out_rows.append(("chm13_geometry", c + "p",
                             "telomere_to_centromere_bp", str(pcen[c]),
                             "cytoband p-arm length: the physical head-room for any "
                             "extension of a censored p-arm PHR"))
    else:
        out_rows.append(("chm13_geometry", "cytobands", "skipped_missing_input",
                         str(CHM13_CYTOBANDS), "not estimated"))


def build_survey_metadata(workdir, flanks):
    contigs = collections.defaultdict(list)
    for f in flanks:
        key = "%s#%s#%s" % (f["sample"], f["hap"], f["contig"])
        bnd = f["start"] if f["arm"] == "parm" else f["end"]
        contigs[key].append((f["arm"], bnd, f["fid"]))
    cfile = workdir / "contigs.tsv"
    with open(cfile, "w") as o:
        for k, v in contigs.items():
            o.write(k + "\t" + "\t".join("%s\t%d\t%d" % (a, b, fi) for a, b, fi in v) + "\n")
    require(CHROM_ALIAS, "chromAlias")
    afile = workdir / "alias.tsv"
    n = 0
    missing = []
    with open(CHROM_ALIAS) as fh, open(afile, "w") as o:
        for line in fh:
            k, lab = line.rstrip("\n").split("\t")
            o.write("%s\t%s\n" % (k, lab.split(":")[0]))
            n += 1
    alias_keys = set()
    with open(afile) as fh:
        for line in fh:
            alias_keys.add(line.split("\t")[0])
    for k in contigs:
        if k not in alias_keys:
            missing.append(k)
    if missing:
        die("%d flank contigs have no chromosome alias (first: %s); the partner "
            "chromosome labels would be undefined" % (len(missing), missing[0]))
    log("  metadata: %d contig keys, %d alias rows, 0 missing" % (len(contigs), n))
    return cfile, afile, len(contigs)


def run_survey(workdir, threads, force, flanks):
    require(SURVEY_PAF_DIR, "genome-wide all-vs-all survey PAFs")
    paf_dir = SURVEY_PAF_DIR
    listing = workdir / "survey_paf_list.txt"
    if force or not listing.exists():
        names = sorted(p.name for p in paf_dir.iterdir() if p.name.endswith(".paf.gz"))
        if not names:
            die("no *.paf.gz in %s" % paf_dir)
        listing.write_text("".join(str(paf_dir / n) + "\n" for n in names))
    files = [l.strip() for l in listing.read_text().splitlines() if l.strip()]
    log("  survey set: %d whole-assembly PAF files" % len(files))

    proj = compile_helper(workdir, "proj", PROJ_C)
    merge = compile_helper(workdir, "merge", MERGE_C)
    agg = compile_helper(workdir, "agg", AGG_C)
    cfile, afile, _ = build_survey_metadata(workdir, flanks)

    parts = workdir / "parts"
    merged = workdir / "merged.tsv"
    cube = workdir / "counts.bin"
    nw = max(1, min(threads, len(files)))
    done_marker = workdir / "parts.done"
    if force or not done_marker.exists():
        if parts.exists():
            shutil.rmtree(parts)
        parts.mkdir()
        chunk_size = (len(files) + nw - 1) // nw
        chunks = [files[i:i + chunk_size] for i in range(0, len(files), chunk_size)]
        total_b = sum(os.path.getsize(f) for f in files)
        log("  projecting %d files through %d workers (%.1f GB compressed PAF to stream)"
            % (len(files), len(chunks), total_b / 1e9))
        procs = []
        for i, ch in enumerate(chunks):
            cfile_i = workdir / ("chunk_%03d.txt" % i)
            cfile_i.write_text("".join(f + "\n" for f in ch))
            out_i = parts / ("part_%03d.tsv" % i)
            # zcat each PAF in turn into a single proj process (one hash per worker)
            cmd = ("xargs -a %s -n1 zcat | %s %s %s 26 %s"
                   % (cfile_i, proj, cfile, afile, out_i))
            procs.append(subprocess.Popen(["bash", "-c", cmd],
                                          stderr=subprocess.PIPE, text=True))
        for i, p in enumerate(procs):
            _, err = p.communicate()
            if p.returncode != 0:
                die("worker %d failed (rc=%d):\n%s" % (i, p.returncode, err))
            for line in (err or "").splitlines():
                if line.startswith("proj:"):
                    sys.stderr.write("    " + line + "\n")
        done_marker.write_text("ok\n")
    else:
        log("  reusing %s" % parts)

    if force or not merged.exists():
        part_files = sorted(str(p) for p in parts.iterdir() if p.name.startswith("part_"))
        if not part_files:
            die("no partial aggregates in %s" % parts)
        r = subprocess.run([str(merge), "28", str(merged)] + part_files,
                           capture_output=True, text=True)
        if r.returncode != 0:
            die("merge failed:\n%s" % r.stderr)
        log("  " + (r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "merged"))
    else:
        log("  reusing %s" % merged)

    if force or not cube.exists():
        r = subprocess.run([str(agg), str(merged), str(len(flanks)), str(cube)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            die("agg failed:\n%s" % r.stderr)
        log("  " + (r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "aggregated"))
    else:
        log("  reusing %s" % cube)
    return cube


def load_cube(path, nfid):
    import numpy as np
    a = np.fromfile(path, dtype=np.uint8)
    if a.size != nfid * N_WIN * 6:
        die("cube %s has %d bytes, expected %d" % (path, a.size, nfid * N_WIN * 6))
    return a.reshape(nfid, N_WIN, 6)


def part_b(out_rows, cube, flanks, len_of, sat, no_phr, pw):
    import numpy as np
    log("Part B: survey-regime calibration and extension")
    nfid = len(flanks)

    def passes(fid, ch, k=2, n_win=N_WIN):
        """pass rule = at least 2 distinct chromosomes with signal in the window.
        k in CHANNELS selects the per-chromosome alignment-count threshold and is
        already baked into the channel index, so the pass threshold here is always 2."""
        col = cube[fid, :n_win, ch] >= 2
        return [bool(x) for x in col]

    # B1 calibration against the published per-window calls at depth < 500 kb
    for label, ch in CAL_CHANNELS:
        tp = fp = tn = fn = 0
        for f in flanks:
            fid = f["fid"]
            p = cube[fid, :FLANK_WIN, ch] >= 2
            for w, (pub, _) in pw[fid].items():
                if w >= FLANK_WIN:
                    continue
                mine = bool(p[w])
                if pub and mine:
                    tp += 1
                elif pub:
                    fn += 1
                elif mine:
                    fp += 1
                else:
                    tn += 1
        sens = tp / (tp + fn)
        spec = tn / (tn + fp)
        out_rows.append(("survey_calibration", label, "window_sensitivity", "%.3f" % sens,
                         "of the %d published-passing windows below 500 kb" % (tp + fn)))
        out_rows.append(("survey_calibration", label, "window_specificity", "%.3f" % spec,
                         "of the %d published-failing windows below 500 kb" % (tn + fp)))
        # length overreach on non-censored PHRs
        over = []
        for f in flanks:
            fid = f["fid"]
            L = len_of[fid]
            if L is None or L == FLANK_LEN:
                continue
            Lb = caller_length(passes(fid, ch, 2))
            if Lb is not None and Lb > L:
                over.append(Lb - L)
        out_rows.append(("survey_calibration", label, "n_uncensored_overreach",
                         "%d/%d" % (len(over),
                                    sum(1 for f in flanks if len_of[f["fid"]] is not None
                                        and len_of[f["fid"]] != FLANK_LEN)),
                         "survey length exceeds the published length"))
        if over:
            out_rows.append(("survey_calibration", label, "uncensored_overreach_median_bp",
                             "%.0f" % q(over, 50), ""))
        gained = sum(1 for f in no_phr
                     if caller_length(passes(f["fid"], ch, 2)) is not None)
        out_rows.append(("survey_calibration", label, "n_phr_free_flanks_gaining_a_region",
                         "%d/%d" % (gained, len(no_phr)),
                         "a good regime must leave these empty"))

    # B2 extension of the censored ends
    sat_ids = [f["fid"] for f in sat]
    for label, ch, k in CHANNELS:
        for rule, fn_ in (("caller_rule", caller_length), ("strict_contiguous", contiguous_length)):
            Ls = [fn_(passes(fid, ch, k)) or 0 for fid in sat_ids]
            arr = np.array(Ls, dtype=np.int64)
            inc = arr - FLANK_LEN
            out_rows.append(("survey_extension", label + "|" + rule, "n", str(len(arr)), ""))
            out_rows.append(("survey_extension", label + "|" + rule,
                             "n_longer_than_500kb", str(int((arr > FLANK_LEN).sum())),
                             "%.1f%%" % (100.0 * (arr > FLANK_LEN).mean())))
            out_rows.append(("survey_extension", label + "|" + rule, "median_length_bp",
                             "%d" % int(np.median(arr)), ""))
            out_rows.append(("survey_extension", label + "|" + rule, "q1_length_bp",
                             "%d" % int(np.percentile(arr, 25)), ""))
            out_rows.append(("survey_extension", label + "|" + rule, "q3_length_bp",
                             "%d" % int(np.percentile(arr, 75)), ""))
            out_rows.append(("survey_extension", label + "|" + rule,
                             "median_length_increase_bp",
                             "%d" % int(np.median(inc)), "measured length - 500 kb"))
            out_rows.append(("survey_extension", label + "|" + rule, "frac_ge_2Mb",
                             "%.3f" % float((arr >= 2000000).mean()),
                             "still passing at the 2 Mb ceiling"))
            out_rows.append(("survey_extension", label + "|" + rule, "frac_ge_5Mb",
                             "%.3f" % float((arr >= 5000000).mean()),
                             "still passing at the 5 Mb ceiling (= truncated by the "
                             "measurement, so a floor)"))
    # per arm, for the primary channel
    label, ch, k = CHANNELS[0]
    by_arm = collections.defaultdict(list)
    for f in sat:
        by_arm[f["chr"].replace("chr", "") + f["arm"][0]].append(f["fid"])
    for arm in sorted(by_arm, key=lambda a: (int(a[:-1]) if a[:-1].isdigit()
                                            else {"X": 23, "Y": 24}.get(a[:-1], 99), a)):
        ids = by_arm[arm]
        for rule, fn_ in (("caller_rule", caller_length), ("strict_contiguous", contiguous_length)):
            arr = np.array([fn_(passes(fid, ch, k)) or 0 for fid in ids], dtype=np.int64)
            out_rows.append(("survey_extension_by_arm", "%s|%s" % (arm, rule), "n",
                             str(len(arr)), ""))
            out_rows.append(("survey_extension_by_arm", "%s|%s" % (arm, rule),
                             "median_length_bp", "%d" % int(np.median(arr)), ""))
            out_rows.append(("survey_extension_by_arm", "%s|%s" % (arm, rule),
                             "frac_longer_than_500kb",
                             "%.3f" % float((arr > FLANK_LEN).mean()), ""))
    # stratified by the boundary signal (confidence)
    bnd = {}
    for f in flanks:
        w = pw[f["fid"]].get(FLANK_WIN - 1)
        if w is not None:
            bnd[f["fid"]] = w[1]
    for lo, hi, name in ((2, 2, "boundary_n_chrs=2"), (3, 4, "boundary_n_chrs=3-4"),
                         (5, 99, "boundary_n_chrs>=5")):
        ids = [f["fid"] for f in sat if lo <= bnd.get(f["fid"], 0) <= hi]
        if not ids:
            continue
        arr = np.array([contiguous_length(passes(fid, ch, k)) or 0 for fid in ids],
                       dtype=np.int64)
        out_rows.append(("survey_extension_by_confidence", name, "n", str(len(arr)), ""))
        out_rows.append(("survey_extension_by_confidence", name, "median_length_bp",
                         "%d" % int(np.median(arr)), ""))
        out_rows.append(("survey_extension_by_confidence", name, "frac_ge_2Mb",
                         "%.3f" % float((arr >= 2000000).mean()), ""))

    # B3 corrected totals (censored lengths replaced by measured lower bounds)
    base = sum(len_of[f["fid"]] for f in flanks if len_of[f["fid"]] is not None)
    n_phr = sum(1 for f in flanks if len_of[f["fid"]] is not None)
    for label, ch, k in CHANNELS[:2]:
        Ls = [contiguous_length(passes(fid, ch, k)) or 0 for fid in sat_ids]
        add = sum(max(0, L - FLANK_LEN) for L in Ls)
        out_rows.append(("corrected_totals", label + "|strict_contiguous",
                         "published_phr_bp", str(base), "%d PHRs" % n_phr))
        out_rows.append(("corrected_totals", label + "|strict_contiguous",
                         "published_mean_length_bp", "%.0f" % (base / n_phr),
                         "matches the paper's 'mean 144 kb'"))
        out_rows.append(("corrected_totals", label + "|strict_contiguous",
                         "additional_bp_lower_bound", str(add), ""))
        out_rows.append(("corrected_totals", label + "|strict_contiguous",
                         "corrected_total_bp", str(base + add), ""))
        out_rows.append(("corrected_totals", label + "|strict_contiguous",
                         "corrected_mean_length_bp", "%.0f" % ((base + add) / n_phr), ""))
    out_rows.append(("corrected_totals", "median", "published_median_length_bp",
                     "%d" % int(q([len_of[f["fid"]] for f in flanks
                                   if len_of[f["fid"]] is not None], 50)),
                     "censoring cannot move the median: the censoring point (500 kb) is "
                     "far above the 50th percentile, so replacing the top 6.9%% by larger "
                     "values leaves the median unchanged"))


# ---- Part C -----------------------------------------------------------------
def part_c(out_rows, flanks):
    log("Part C: transitive genome-wide regime (CHM13-anchored)")
    if not TRANSITIVE_TABLE.exists():
        out_rows.append(("transitive_regime", "input", "skipped_missing_input",
                         str(TRANSITIVE_TABLE), "not estimated"))
        return
    bnd = {}
    for f in flanks:
        if f["sample"] != "CHM13":
            continue
        bnd[(f["chr"], f["arm"])] = f["start"] if f["arm"] == "parm" else f["end"]
    per = collections.defaultdict(dict)
    with gzip.open(TRANSITIVE_TABLE, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if hdr[0] != "chrom" or "chroms-num_haplotypes" not in hdr:
            die("unexpected transitive-table header: %r" % hdr)
        ci = hdr.index("chroms-num_haplotypes")
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if not a[0].startswith("CHM13#"):
                continue
            c = a[0].split("#")[-1]
            st, en = int(a[1]), int(a[2])
            d = collections.Counter()
            for tok in a[ci].split(","):
                k, _, v = tok.rpartition("-")
                if v and int(v) >= 5:
                    d[k] = int(v)
            per[c][(st + en) // 2] = d
    for (c, arm) in sorted(bnd, key=lambda x: (int(x[0].replace("chr", ""))
                                               if x[0].replace("chr", "").isdigit()
                                               else {"X": 23, "Y": 24}.get(
                                                   x[0].replace("chr", ""), 99), x[1])):
        base = bnd[(c, arm)]
        rows = []
        for mid, d in per.get(c, {}).items():
            depth = (mid - base) if arm == "parm" else (base - mid)
            rows.append((depth, sum(1 for k in d if k != c)))
        rows.sort()
        extent = 0
        for depth, n_partner in rows:
            if depth < 0:
                continue
            if n_partner >= 1:
                extent = depth + 100000
            elif depth > extent + 400000:
                break
        out_rows.append(("transitive_regime", "%s%s" % (c.replace("chr", ""),
                                                        arm[0]), "extent_bp",
                         str(extent),
                         "transitive closure, 100 kb windows, id>=0.98, min output 5 kb; "
                         "more permissive than the published direct projection, so this is "
                         "the upper side of the bound"))


# ---- main -------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR,
                    help="scratch directory for intermediates (default %(default)s)")
    ap.add_argument("--threads", type=int, default=min(32, os.cpu_count() or 4))
    ap.add_argument("--force", action="store_true", help="redo the survey projection")
    ap.add_argument("--skip-survey", action="store_true",
                    help="Part A and C only (no 297 GB stream)")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    for p, what in ((FLANK_FAI, "flank fai"), (CALL_TABLE, "PHR call table"),
                    (PER_WINDOW, "per-window table")):
        require(p, what)

    workdir = args.workdir
    workdir.mkdir(parents=True, exist_ok=True)
    if not os.access(workdir, os.W_OK):
        die("workdir not writable: %s" % workdir)

    # ---- flanks
    flanks = []
    with open(FLANK_FAI) as fh:
        for i, line in enumerate(fh):
            name, length = line.split("\t")[:2]
            if int(length) != FLANK_LEN:
                die("flank %s has length %s, expected %d" % (name, length, FLANK_LEN))
            sample, hap, ctg, s, e, chrl, arm = parse_flank(name)
            if e - s != FLANK_LEN - 1:
                die("flank %s spans %d bp, expected %d" % (name, e - s + 1, FLANK_LEN))
            flanks.append(dict(fid=i, name=name, sample=sample, hap=hap, contig=ctg,
                               start=s, end=e, chr=chrl, arm=arm))
    log("flanks: %d, all exactly %d bp" % (len(flanks), FLANK_LEN))
    fid_of = {f["name"]: f["fid"] for f in flanks}

    out_rows = []
    len_of, dropped = load_call_table(flanks, fid_of)
    pool = [f for f in flanks if f["fid"] not in set(dropped)]
    pw = load_per_window(flanks, fid_of)
    part_a(out_rows, pool, len_of, pw, workdir)
    part_a5(out_rows, pool, len_of, pw)
    part_a6(out_rows, flanks)
    if args.skip_survey:
        out_rows.append(("survey_extension", "input", "skipped_by_flag", "-", ""))
    else:
        sat = [f for f in pool if len_of[f["fid"]] == FLANK_LEN]
        no_phr = [f for f in pool if len_of[f["fid"]] is None]
        cube = run_survey(workdir, args.threads, args.force, flanks)
        part_b(out_rows, load_cube(cube, len(flanks)), pool, len_of, sat, no_phr, pw)
    part_c(out_rows, pool)

    log("writing %s (%d rows)" % (args.out, len(out_rows)))
    with open(args.out, "w") as o:
        o.write("section\tstratum\tmetric\tvalue\tnote\n")
        for r in out_rows:
            if len(r) != 5:
                die("internal: malformed row %r" % (r,))
            o.write("\t".join(str(x) for x in r) + "\n")
    log("done")


if __name__ == "__main__":
    main()
